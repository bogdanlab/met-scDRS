# met_scdrs/cli.py

import argparse
from met_scdrs.core import score_cells

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required = True)
    args = parser.parse_args()

    print(f"[Computing]")
    score_cells(args.input)