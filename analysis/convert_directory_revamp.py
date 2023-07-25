import numpy as np
import uproot as ur
from typing import Optional,Union,List,Dict

from os import listdir, mkdir, getcwd
from os.path import isfile, join, basename, isdir

# TODO: Support slicing together of multiple ROOT files based on matching columns (e.g. number of run, event)

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

# Loads a LCIO file and converts it into a numpy representation
def lcio_to_numpy(source_file: str, event_nr: int, collection: str):
    ...

# Converts the first file in source_files and adds columns of all subsequent files by joining join_by columns and merging according to columns_to_merge
# source_files     : list of absolute paths
# in_file_location : name of TTree (or list, matching per source_file)
# join_by          : columns by which data of different files is joined
# columns_to_merge : for each entry in source_files, a dict of form source_column: target_column, ...
def root_multiple_to_numpy(source_files: List[str],
                  in_file_location: Union[str,List[str]],
                  join_by: List[str],
                  columns_to_merge: List[Dict[str, str]]):
    
    # Import complete first file
    npa = root_to_numpy(source_files[0], in_file_location[0] if isinstance(in_file_location, list) else in_file_location)
    
    # Get correct column names and types for conversion
    dtype_names = list(npa.dtypes.keys())
    dtype_join_by = []
    
    for join_by_column in join_by:
        dtype_join_by.append((join_by_column, npa.dtypes[join_by_column].name))
    
    for j in range(1, len(source_files)):
        in_file_location = in_file_location[j] if isinstance(in_file_location, list) else in_file_location
        source_path = source_files[j]
        merge_columns = columns_to_merge[j]
        
        source_columns = []
        target_columns = []
        
        for source_column in merge_columns:
            target_column = merge_columns[source_column]
            if target_column in dtype_names:
                raise Exception("Error: {0} already exists in target!".format(target_column))
            
            if source_column in source_columns:
                raise Exception("Error: {0} already exists in source!".format(source_column))
            
            source_columns.append(source_column)
            target_columns.append(target_column)

        with ur.open(source_path) as file:
            # Find data in file
            data = file[in_file_location]
            dtype_names = data.typenames()
            
            # Fetch join_by columns in current file
            join_by_columns = np.zeros(data.num_entries, dtype=dtype_join_by)
            for join_by_col in join_by:
                join_by_columns[join_by_col] = data[join_by_col].array()

            for i in range(0, len(source_columns)):
                source_column = source_columns[i]
                
                col_values_src = np.array(data[source_column].array(), dtype=convert_type(dtype_names[source_column]))
                col_values_tgt = np.zeros(len(npa)                   , dtype=convert_type(dtype_names[source_column]))
            
                join_by_a = npa[join_by]
                join_by_b = join_by_columns
                
                join_by_a_view = join_by_a.view([('',join_by_a.dtype)]*len(join_by_a.dtype.names))
                join_by_b_view = join_by_b.view([('',join_by_b.dtype)]*len(join_by_b.dtype.names))
                
                intersection, a_idx, b_idx = np.intersect1d(join_by_a_view, join_by_b_view, return_indices=True)
                
                for i in range(0, len(a_idx)):
                    col_values_tgt[a_idx[i]] = col_values_src[b_idx[i]]
                    
                npa = np.append(npa, col_values_tgt, axis=1)

    return npa

# Loads a ROOT file and converts it into a numpy representation
def root_to_numpy(source_path: str,
                  in_file_location: str,
                  merge_with_np_array: Optional[np.ndarray] = None,
                  join_by: Optional[list]=None,
                  merge_columns: Optional[list]=None) -> np.ndarray:
    
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
            print("Skipping file <" + bname + "> (exists)")
        else:
            convert_file(join(source_path, filename), in_file_location, output_path)
            n_converted += 1
            
    return n_converted

#Test
root_multiple_to_numpy(["/nfs/dust/ilc/user/bliewert/fullflow_v3/comparison/root/prod/compare_reco.root",
                        "/nfs/dust/ilc/user/bliewert/fullflow_v3/zhh/root/rv02-02-03.sv02-02-03.mILD_l5_o1_v02_nobg.E500-TDR_ws.I403001.Pe2e2hh.eL.pR.n000.d_dstm_15807_0_MisclusteringSig.root"],
                       ["dataTree", "eventTree"],
                       ["run", "event"],
                       [
                           {
                               "asd": ""
                           }
                       ])