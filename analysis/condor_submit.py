# 2024, Bryan Bliewert
# bryan.bliewert@desy.de

# Example usage: python analysis/condor_submit.py jobs/test.sub -n 4

import argparse
import os
import subprocess
import threading
import re
from tqdm.auto import tqdm
from queue import Queue
from math import floor
from typing import Optional, Callable

parser = argparse.ArgumentParser(
                    prog='condor_submit',
                    description='A shim to run condor submit files on a local system with n threads. Supported fields: executable, arguments, queue from file, log, output and error. Supported arguments: $(file) $(ProcId).',
                    epilog='Not for use on batch farms')

parser.add_argument('filename', help='Full path to the .sub file')
parser.add_argument('-n', '--nthreads', help='Number of threads to spawn', default=4)

args = parser.parse_args()
print(args)

# Lowercase values to identify
ARG_IDS = {
    'FILE': ['$(file)'],
    'PROCID': ['$(ProcId)', '$(Process)']
}

# From https://stackoverflow.com/questions/1717393/is-the-operator-thread-safe-in-python
class ThreadSafeCounter():
    def __init__(self):
        self.lock = threading.Lock()
        self.value = 0

    def increment(self):
        with self.lock:
            self.value += 1

    def decrement(self):
        with self.lock:
            self.value -= 1
            
class ThreadSafePBar():
    def __init__(self, total:int, desc:str='Total progress', position:int=0):
        self.lock = threading.Lock()
        self.pbar = tqdm(total=total, desc=desc, position=position)

    def update(self):
        with self.lock:
            self.pbar.update()

# From https://stackoverflow.com/questions/2581817/python-subprocess-callback-when-cmd-exits
def popen_and_call(on_exit:Callable, popen_args):
    """
    Runs the given args in a subprocess.Popen, and then calls the function
    on_exit when the subprocess completes.
    on_exit is a callable object, and popen_args is a list/tuple of args that 
    would give to subprocess.Popen.
    """
    def run_in_thread(on_exit, popen_args):
        proc = subprocess.Popen(**popen_args)
        proc.wait()
        on_exit()
        return
    
    thread = threading.Thread(target=run_in_thread, args=(on_exit, popen_args))
    thread.start()
    # returns immediately after the thread starts
    return thread

def start_queue(executable:str, arguments:list, files:list, streams:dict, nthreads=args.nthreads):
    file_queue = Queue()
    for file in files:
        file_queue.put(file)
    
    scope = {
        'counter': ThreadSafeCounter(),
        'n_total': len(files),
        'file_queue': file_queue,
        'progress': ThreadSafePBar(total=len(files))
    }
    
    # Add progress bars
    nthreads_spawn = min(int(nthreads), scope['n_total'])
    
    queue_sizes = []
    even = floor(len(files) / nthreads_spawn)
    remainder = len(files) % nthreads_spawn
    
    for i in range(nthreads_spawn):
        queue_size = even + (1 if remainder > i else 0)
        queue_sizes.append(queue_sizes)
        scope[f'progress_{str(i)}'] = ThreadSafePBar(total=queue_size, position=i, desc=f'Thread [{str(i)}]')
    
    # Define callbacks    
    def finalize():
        pass
        #if scope['counter'].value == scope['n_total']:
        #    print('Finished execution')
    
    def on_done(executable, arguments, streams_in:dict, scope:dict, thread_id:int):
        def result():
            scope[f'progress_{str(thread_id)}'].update()
            next(executable, arguments, streams_in, scope, thread_id)
            
        return result
    
    def next(executable, arguments, streams_in:dict, scope:dict, thread_id:int):
        if scope['file_queue'].empty():
            return finalize()
        
        file = scope['file_queue'].get()
        scope['counter'].increment()
        procid = scope['counter'].value - 1

        # Prepare streams (if requires)
        streams = { **{ 'error': '/dev/null', 'output': '/dev/null' }, **streams_in }
        for stream_key, stream in streams_in.items():
            for key in ARG_IDS:
                choices = ARG_IDS[key]
                val = str(procid) if key == 'PROCID' else file
                
                for choice in choices:
                    choice_ci = re.compile(re.escape(choice), re.IGNORECASE)
                    stream = choice_ci.sub(val, stream)
            
            streams[stream_key] = stream
        
        # Compile arguments
        compiled_arguments = []
        for arg_key in arguments:
            if arg_key == 'FILE':
                compiled_arguments.append(file)
            elif arg_key == 'PROCID':
                compiled_arguments.append(str(procid))
        
        #executable = "sleep"
        #compiled_arguments = ["5"]
        
        if True:
            stdout = open(streams['output'], 'w')
            stderr = open(streams['error'], 'w')
            
            popen_and_call(on_done(executable, arguments, streams_in, scope, thread_id),
                           { 'args': [executable, *compiled_arguments],  'stdout': stdout, 'stderr': stderr })
    
    # Initial start
    for i in range(nthreads_spawn):        
        next(executable, arguments.copy(), streams.copy(), scope, i)
    
    #print(executable, arguments)

def parse_submit_file(args):
    # Options to be populated
    executable = ''
    arguments = []
    files = []
    streams = {}
    with open(args.filename) as file:    
        for line in file:
            # Remove whitespaces around line
            line = line.strip()
            
            if line.startswith('executable'):
                executable_candidate = line.split('=')[1].strip()
                if os.path.isfile(executable_candidate):
                    executable = executable_candidate
                else:
                    raise Exception(f'Invalid executable <{executable_candidate}>')
            
            if line.startswith('arguments'):
                arguments_string = line.split('=')[1].strip()
                for arg in arguments_string.split(' '):
                    for key in ARG_IDS:
                        choices = ARG_IDS[key]
                        for choice in choices:
                            if arg.lower() == choice.lower():
                                arguments.append(key)
                
                if not len(arguments):
                    raise Exception(f'Could not read arguments from <{arguments_string}>')
            
            if line.startswith('queue file from '):
                queue_path = line.split('queue file from ')[1].strip()
                if os.path.isfile(queue_path):
                    with open(queue_path) as queue_file:
                        for qline in queue_file:
                            qline = qline.strip()
                            files.append(qline)
                else:
                    raise Exception(f'Could not read queue from file <{queue_path}>')
            
            for key in ['log', 'error', 'output']:
                if line.startswith(key):
                    stream_candidate = line.split('=')[1].strip()
                    stream_candidate_dir = os.path.dirname(stream_candidate)
                    if os.path.isdir(stream_candidate_dir):
                        streams[key] = stream_candidate
                    else:
                        raise Exception(f'Error: Parent directory <{stream_candidate_dir}> does not exist')
    
    if executable == '' or not len(arguments) or not len(files):
        raise Exception('Invalid job specification')
    
    return executable, arguments, files, streams

start_queue(*parse_submit_file(args=args))