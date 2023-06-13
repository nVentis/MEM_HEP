from setuptools import setup

setup(
    name='mem_hep_ana',
    version='0.1.0',
    py_modules=['convert_directory'],
    install_requires=[
        'Click',
    ],
    entry_points={
        'console_scripts': [
            'convert_directory = convert_directory:cli',
        ],
    },
)