def module_reload(prefix:str):
    import importlib
    import sys

    for k,v in list(sys.modules.items()):
        if k.startswith(prefix):
            importlib.reload(v)