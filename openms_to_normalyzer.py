#!/usr/bin/env python3

import argparse
import pandas as pd

JOINT_COLS = ['rt', 'mz', 'intensity', 'charge', 'quality']
SAMPLE_COLS = ['intensity']

VERSION = '1.0.0'


def main():

    args = parse_arguments()

    cons_lines = get_consensus_lines(args.input)
    cons_fp = '{}.cons_only'.format(args.input)
    generate_consensus_intermediate(cons_fp, cons_lines)

    ms_df = pd.read_csv(cons_fp, sep=args.delim_in)
    print(ms_df.head())

    design_matrix = pd.read_csv(args.design, sep="\s+")

    normalyzer_df = setup_normalyzer_df(ms_df, design_matrix)

    print("Writing dataframe with shape {}, to {}"
          .format(normalyzer_df.shape, args.output))
    normalyzer_df.to_csv(args.output, sep=args.delim_out, header=False, index=False, na_rep='NA')
    print("Done!")


def get_consensus_lines(in_fp):

    """Extract consensus lines from raw OpenMS file"""

    cons_lines = list()
    with open(in_fp) as in_fh:
        for line in in_fh:
            line = line.rstrip()
            if line.startswith('CONSENSUS') or line.startswith('#CONSENSUS'):
                cons_lines.append(line)
    return cons_lines


def generate_consensus_intermediate(cons_fp, cons_lines):

    """Extract consensus rows from OpenMS file and write to new file"""

    print('TODO: Omit intermediate file')

    with open(cons_fp, 'w') as out_fh:
        for line in cons_lines:
            print(line, file=out_fh)


def get_sample_numbers(cons_df):

    """Extract sample numbers from OpenMS header"""

    all_cols = list(cons_df.columns)
    sample_nbrs = set([int(col.split('_')[-1]) for col in all_cols[1:] if not col.endswith('cf')])
    return sample_nbrs


def setup_normalyzer_df(ms_df, design_matrix, nan_fill='NA'):

    """
    Setup Normalyzer dataframe. The dataframe consists of two-row header, and databody
    It consists of both annotation columns, and sample columns
    """

    target_sample_cols = get_sample_colnames(ms_df)
    target_annot_cols = get_annot_colnames()

    normalyzer_vals = pd.DataFrame(ms_df[target_annot_cols + target_sample_cols].fillna(nan_fill).applymap(str))

    headers = setup_normalyzer_header(design_matrix, target_annot_cols, normalyzer_vals)
    normalyzer_df = pd.concat([headers, normalyzer_vals])

    return normalyzer_df


def get_sample_colnames(ms_df):

    """
    Generate names for sample columns
    """

    sample_numbers = get_sample_numbers(ms_df)

    target_sample_cols = list()
    for sample in sample_numbers:
        for col in SAMPLE_COLS:
            target_sample_cols.append('{attr}_{sample}'.format(attr=col, sample=sample))
    return target_sample_cols


def get_annot_colnames():

    """
    Retrieve names for OpenMS columns
    """

    target_annot_cols = list()
    for col in JOINT_COLS:
        target_annot_cols.append('{}_{}'.format(col, 'cf'))
    return target_annot_cols


def setup_normalyzer_header(design_matrix, target_annot_cols, normalyzer_vals):

    """
    Setup two top rows in Normalyzer matrix
    """

    sample_head = [-1] + [0] * (len(target_annot_cols)-1) + list(design_matrix['biorep'])
    sample_head_str = [str(e) for e in sample_head]

    nbr_annot_cols = len(target_annot_cols)
    label_row = list(normalyzer_vals.columns)[:nbr_annot_cols] + list(design_matrix['name'])

    headers = pd.DataFrame([sample_head_str, label_row])
    headers.columns = normalyzer_vals.columns

    return headers


def parse_arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', help='OpenMS TSV report', required=True)
    parser.add_argument('-o', '--output', help='Normalyzer formatted report', required=True)
    parser.add_argument('--design', help='Sample design file', required=True)

    parser.add_argument('--delim_in', default='\t')
    parser.add_argument('--delim_out', default='\t')
    parser.add_argument('-v', '--version', action='version', version=VERSION)

    return parser.parse_args()


if __name__ == '__main__':
    main()
