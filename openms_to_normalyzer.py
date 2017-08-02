#!/usr/bin/env python3

from typing import List, Set

import argparse
import pandas as pd
import numpy as np

DF = pd.DataFrame

JOINT_COLS = ['#rt', 'mz', 'intensity', 'charge', 'quality']
SAMPLE_COLS = ['intensity']
ANNOT_COLS = ['peptide_0', 'protein_0']
VERSION = '2.0.0'


def main():

    args = parse_arguments()

    ms_df = pd.read_csv(args.input, sep=args.delim_in, skiprows=2)
    design_matrix = pd.read_csv(args.design, sep="\s+")
    normalyzer_df = setup_normalyzer_df(ms_df, design_matrix, require_annot=args.require_annot)

    # if args.require_annot:
    #     normalyzer_df = normalyzer_df.loc[normalyzer_df['peptide_0'] != np.nan]

    print("Writing dataframe with shape {}, to {}"
          .format(normalyzer_df.shape, args.output))
    normalyzer_df.to_csv(args.output, sep=args.delim_out, header=False, index=False, na_rep='NA')
    print("Done!")


def setup_normalyzer_df(ms_df: DF, design_matrix: DF, require_annot=False, nan_fill='NA') -> DF:

    """
    Setup Normalyzer dataframe. The dataframe consists of two-row header, and databody
    It consists of both annotation columns, and sample columns
    """

    target_sample_cols = get_sample_colnames(ms_df)
    target_annot_cols = get_annot_colnames()

    normalyzer_vals = pd.DataFrame(ms_df[target_annot_cols + target_sample_cols].fillna(nan_fill).applymap(str))

    if require_annot:
        normalyzer_vals = normalyzer_vals.loc[normalyzer_vals['peptide_0'] != nan_fill]

    headers = setup_normalyzer_header(design_matrix, target_annot_cols, normalyzer_vals)
    normalyzer_df = headers.append(normalyzer_vals)

    return normalyzer_df


def get_sample_colnames(ms_df: DF) -> List[str]:

    """
    Generate names for sample columns
    """

    sample_numbers = get_sample_numbers(ms_df)

    target_sample_cols = list()
    for sample in sample_numbers:
        for col in SAMPLE_COLS:
            target_sample_cols.append('{attr}_{sample}'.format(attr=col, sample=sample))
    return target_sample_cols


def get_sample_numbers(cons_df: pd.DataFrame) -> Set[int]:

    """Extract sample numbers from OpenMS header"""

    all_cols = list(cons_df.columns)
    sample_nbrs = set([int(col.split('_')[-1]) for col in all_cols[1:] if not col.endswith('cf')])
    return sample_nbrs


def get_annot_colnames() -> List[str]:

    """
    Retrieve names for OpenMS columns
    """

    target_annot_cols = list()

    for col in JOINT_COLS:
        target_annot_cols.append('{}_{}'.format(col, 'cf'))

    target_annot_cols += ANNOT_COLS

    return target_annot_cols


def setup_normalyzer_header(design_matrix: DF, annot_cols: List[str], normalyzer_vals:DF) -> DF:

    """
    Setup two top rows in Normalyzer matrix
    """

    # Get numbers set up as list of stringified numbers ('-1', '0', '0', '1', '1')
    nbr_annot_cols = len(annot_cols)
    sample_head = [-1] + [0] * (nbr_annot_cols - 1) + list(design_matrix['biorepgroup'])
    sample_head_str = [str(e) for e in sample_head]

    # Get text-information about each column
    label_row = list(normalyzer_vals.columns)[:nbr_annot_cols] + list(design_matrix['name'])

    headers = pd.DataFrame([sample_head_str, label_row])
    headers.columns = normalyzer_vals.columns

    return headers


def parse_arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', help='OpenMS TSV report', required=True)
    parser.add_argument('-o', '--output', help='Normalyzer formatted report', required=True)
    parser.add_argument('--design', help='Sample design file', required=True)

    parser.add_argument('--require_annot', help='Output only features with matched peptide', action='store_true',
                        default=False)

    parser.add_argument('--delim_in', default='\t')
    parser.add_argument('--delim_out', default='\t')
    parser.add_argument('-v', '--version', action='version', version=VERSION)

    return parser.parse_args()


if __name__ == '__main__':
    main()
