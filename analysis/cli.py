import sys
import argparse

# Available packages
import convert_directory

if len(sys.argv) > 1 & len(sys.argv[1]) > 0:
    print(f"Arguments count: {len(sys.argv)}")
    for i, arg in enumerate(sys.argv):
        print(f"Argument {i:>6}: {arg}")

