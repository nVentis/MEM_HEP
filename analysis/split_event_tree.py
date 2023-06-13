import numpy as np
import pandas as pd

# Splits a dataframe with separate columns for reco and true properties to one with the properties and a different label "type" with values "true" and "reco" 
def split_event_tree(df,
                     properties=["sigma", "nll"],
                     props_shared=["run", "event", "is_zhh", "is_zzh", "h1z2_decay_pdg", "llr"],
                     type_names=["zhh", "zzh"],
                     ttype_cols=["is_zhh", "is_zzh"], # true type columns
                     column_concat: str ="_"):
    #type_dict = {}
    #for property in properties:
    #    type_dict[property] = pd.Series(dtype=df.dtypes[type_names[0] + column_concat + property])

    # Rename
    rename_dict = {}
    for property in properties:
        rename_dict[type_names[0] + column_concat + property] = property

    df = df.rename(columns=rename_dict)
    df["mem_type"] = type_names[0]

    # Slice and re-add
    for type in type_names:
        if type == type_names[0]:
            continue

        #column_names = [word.replace(type + column_concat, "") for word in ["asd_123", "asd_345"]]
        props_type_specific = list(map(lambda cn: type+column_concat+cn, properties))
        props_slice = props_shared + props_type_specific
        sliced = df.loc[:, props_slice]

        df.drop(columns=props_type_specific, inplace=True)

        rename_dict = {}
        for property in properties:
            rename_dict[type + column_concat + property] = property

        sliced.rename(columns=rename_dict, inplace=True)
        sliced["mem_type"] = type        

        df = pd.concat([df, sliced])

    # Apply true_type column
    condlist = []
    for col in ttype_cols:
        condlist.append(df[col].astype(bool))

    df["true_type"] = np.select(condlist, type_names, default="none")

    df.reset_index(drop=True, inplace=True)

    return df

def ttype_column(df, column="true_type", type_names = ["zhh", "zzh"], check = ["is_zhh", "is_zzh"]):
    condlist = []
    for col in check:
        condlist.append(df[col].astype(bool))

    df[column] = np.select(condlist, type_names, default="none")
    df.reset_index(drop=True, inplace=True)

    return df