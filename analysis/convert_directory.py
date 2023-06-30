#import click
import numpy as np
import uproot as ur
from typing import Optional,Union

from os import listdir, mkdir, getcwd
from os.path import isfile, join, basename, isdir

type_map = {
    "int32_t": "i4",
    "float": "float32",
    "double": "float64"
}

def convert_type(type_name: str):
    if type_name in type_map:
        return type_map[type_name]
    else:
        return type_name

# Loads a ROOT file and converts it into a numpy representation
def root_to_numpy(source_path: str,
                  in_file_location: str,
                  merge_with_np_array: Optional[np.ndarray] = None,
                  join_by: Optional[list]=None,
                  merge_columns: Optional[list]=None):
    
    with ur.open(source_path) as file:
        # Find data in file
        data = file[in_file_location]
        keys = data.keys()

        # Get correct column names and types for conversion
        dtype_arr = []
        dtype_names = data.typenames()

        for field_name in dtype_names:
            dtype_arr.append((field_name, convert_type(dtype_names[field_name])))

        dtype_list2 = []
        if merge_with_np_array is not None and join_by is not None:
            dtype_list = list(dtype_names.keys())
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

    return out

def convert_file(source_path: str, in_file_location: str, output_path: str):    
    np.save(output_path, root_to_numpy(source_path, in_file_location), allow_pickle=True)

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