# Bryan Bliewert, DESY/TUM
# 07.07.2023

import numpy as np
import uproot as ur
from typing import Optional,Union,List,Dict

from os import listdir, mkdir, getcwd, remove, makedirs
from os.path import isfile, join, basename, isdir, splitext, dirname

type_map = {
    "int32_t": "i4",
    "float": "float32",
    "double": "float64",
    "std::vector<float>": ""
}

def convert_type(type_name: str):
    # Filter out non-supported data-types
    if type_name.startswith("std::vector"):
        return None
    
    if type_name in type_map:
        return type_map[type_name]
    else:
        return type_name

# Loads a LCIO file and converts it into a numpy representation
def lcio_to_numpy(source_file: str, event_nr: int, collection: str):
    ...

# Loads a ROOT file and converts it into a numpy representation
def root_to_numpy(source_path: str,
                  in_file_location: str,
                  merge_with_np_array: Optional[np.ndarray] = None,
                  join_by: Optional[list]=None,
                  merge_columns: Optional[list]=None,
                  null_on_not_found:bool=False) -> Union[np.ndarray, None]:
    
    out = None
    
    with ur.open(source_path) as file:
        # Find data in file
        if in_file_location in file:
            data = file[in_file_location]
            keys = data.keys()

            # Get correct column names and types for conversion
            dtype_arr = []
            dtype_names = data.typenames()
            
            dtype_names_accepted = []

            for key in dtype_names:
                conv_type = convert_type(dtype_names[key])
                if conv_type is not None:
                    dtype_names_accepted.append(key)
                    dtype_arr.append((key, conv_type))

            dtype_list2 = []
            if merge_with_np_array is not None and join_by is not None:
                dtype_list = list(dtype_names_accepted)
                # Get keys only in merge_with_np_array
                dtype_list2 = list(set((merge_columns if merge_columns is not None else list(merge_with_np_array.dtype.fields.keys()))) - set(dtype_list))
                
                for field_name in dtype_list2:
                    dtype_arr.append((field_name, merge_with_np_array.dtype[field_name].name))

            # Convert data to (column-wise) arrays using numpy
            out = np.zeros(data.num_entries, dtype=dtype_arr)

            for i in range(0, len(keys)):
                key = keys[i]
                out[key] = data[key].array()
            
            if merge_with_np_array is not None and join_by is not None:
                join_by_a = out[join_by]
                join_by_b = merge_with_np_array[join_by]
                
                join_by_a_view = join_by_a.view([('',join_by_a.dtype)]*len(join_by_a.dtype.names))
                join_by_b_view = join_by_b.view([('',join_by_b.dtype)]*len(join_by_b.dtype.names))
                
                intersection, a_idx, b_idx = np.intersect1d(join_by_a_view, join_by_b_view, return_indices=True)
                
                for i in range(0, len(a_idx)):
                    out[dtype_list2][a_idx[i]] = tuple(merge_with_np_array[dtype_list2][b_idx[i]])
    
    if out is None:
        if null_on_not_found:
            return out
        else:
            raise Exception(f'{in_file_location} does not exist in file')

    return out

def convert_file(source_path: str, in_file_location: str, output_path:Optional[str] = None)->str:
    """Converts a ROOT to NPY files

    Args:
        source_path (str): file path
        in_file_location (str): name of TTree
        output_path (str): file path or directory; if None, will be equal to source_path

    Returns:
        str: output file path
    """
    
    output_path = dirname(source_path) if output_path is None else output_path
    root_file = basename(source_path)
    
    cnv_file = output_path if splitext(output_path)[1] == ".npy" else join(output_path, splitext(root_file)[0] + ".npy")

    dst_dir = dirname(cnv_file)
    
    if not isdir(dst_dir):
        makedirs(dst_dir)
    
    for file in [cnv_file]:
        if isfile(file):
            remove(file)
            
    np.save(output_path, root_to_numpy(source_path, in_file_location), allow_pickle=True)
    
    return cnv_file

# See https://uproot.readthedocs.io/en/latest/uproot.behaviors.TBranch.iterate.html
def convert_directory(source_path: str, in_file_location: str, output_dir_abs: str = "", overwrite=False):
    """Converts all ROOT trees at position in_file_location of the .root files contained in source_path to npy format in output_dir_abs"""
    dir_contents = listdir(source_path)
    root_files = filter(lambda filename: filename.endswith(".root"), dir_contents)
    
    n_converted = 0

    output_dir = source_path
    if output_dir_abs != "":
        output_dir = output_dir_abs

    if not isdir(output_dir_abs):
        mkdir(output_dir_abs)

    for filename in root_files:
        bname = basename(filename)
        output_path = join(output_dir, bname + ".npy")
        if isfile(output_path):
            if overwrite:
                remove(output_path)
            else:
                print("Skipping file <" + bname + "> (exists)")
        else:
            convert_file(join(source_path, filename), in_file_location, output_path)
            n_converted += 1
            
    return n_converted