import click

@click.group(name="plot")
def cli_plot():
    """Plotting related commands"""
    pass

if __name__ == '__main__':
    cli_plot()