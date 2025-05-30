# met_scdrs/cli.py

import argparse
import met_scdrs

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required = True)
    args = parser.parse_args()

    print(f"[Computing]")
    met_scdrs.score_cells(args.input)