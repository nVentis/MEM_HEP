import click
from analysis import logger

@click.group(name="convert")
def cli_convert():
    """Commands for conversion of specialized ROOT files to numpy"""
    pass

@cli_convert.command()
@click.argument("src")
@click.argument("dst")
@click.option("--tree", default="dataTree", help="Name of TTree in ROOT file")
def file_command(src:str, dst:str, tree:str):
    """Converts the .root file at src to .npy at dst, latter can be file or directory. Skips unsupported types (any vectors)"""
    from analysis.convert import convert_file
                
    convert_file(src, tree, dst)
    logger.info("Converted file {} to {}".format(src, dst))


@cli_convert.command()
@click.argument("src")
@click.argument("dst")
@click.option("--tree", default="dataTree", help="Name of TTree in ROOT file")
@click.option("--overwrite", is_flag=True, default=False, help="Whether or not to overwrite existing files")
def directory_command(src:str, dst:str, tree:str, overwrite:bool):
    """Converts all .root files at src to .npy in directory dst"""
    from analysis.convert import convert_directory
    
    convert_directory(src, tree, dst, overwrite)
    logger.info("Converted directory {} to {}".format(src, dst))
    

if __name__ == '__main__':
    cli_convert()