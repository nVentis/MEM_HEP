#import click
import numpy as np
import uproot as ur

from os import listdir, mkdir, getcwd
from os.path import isfile, join, basename, isdir

type_map = {
    "int32_t": "i4"
}

def convert_type(type_name: str):
    if type_name in type_map:
        return type_map[type_name]
    else:
        return type_name

# Loads a ROOT file and converts it into a numpy representation
def root_to_numpy(source_path, in_file_location, merge_with_np_array=None, merge_by_columns=None):
    with ur.open(source_path) as file:
        # Find data in file
        data = file[in_file_location]
        keys = data.keys()

        # Get correct column names and types for conversion
        dtype_arr = []
        dtype_names = data.typenames()

        for field_name in dtype_names:
            dtype_arr.append((field_name, convert_type(dtype_names[field_name])))
            
        if merge_with_np_array is not None and merge_by_columns is not None:
            dtype_names2 = merge_with_np_array.typenames()
            keys = keys + list(set(dtype_names2) - set(keys))
            
            for field_name in dtype_names2:
                if field_name not in dtype_names:
                    dtype_arr.append((field_name, convert_type(dtype_names[field_name])))

        # Convert data to (column-wise) arrays using numpy
        out = np.zeros(data.num_entries, dtype=dtype_arr)
        out[:] = np.nan

        for i in range(0, len(keys)):
            key = keys[i]
            out[key] = data[key].array()
            
        if merge_with_np_array is not None and merge_by_columns is not None:
            ...

    return out

# merge_with can be a dataframe or numpy array that will be appended to the output
def convert_file(source_path, in_file_location, output_path = "", merge_with=None):
    out = root_to_numpy(source_path, in_file_location, output_path)
    
    if output_path != "":
        np.save(output_path, out, allow_pickle=True)

    return out

# See https://uproot.readthedocs.io/en/latest/uproot.behaviors.TBranch.iterate.html
#@click.command()
#@click.argument("source_path")
#@click.argument("in_file_location")
#@click.option("--output_dir_abs", default="", help="Output directory of the converted .npy files. Per default, the current working directory")
def convert_directory(source_path: str, in_file_location: str, output_dir_abs: str = ""):
    """Converts all ROOT trees at position in_file_location of the .root files contained in source_path to npy format in output_dir_abs"""
    dir_contents = listdir(source_path)
    root_files = filter(lambda filename: filename.endswith(".root"), dir_contents)

    output_dir = source_path
    if output_dir_abs != "":
        output_dir = output_dir_abs

    if not isdir(output_dir_abs):
        mkdir(output_dir_abs)

    for filename in root_files:
        bname = basename(filename)
        output_path = join(output_dir, bname + ".npy")
        if isfile(output_path):
            print("Skipping file <" + bname + "> (exists)")
        else:
            convert_file(join(source_path, filename), in_file_location, output_path)