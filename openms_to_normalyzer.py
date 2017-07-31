#!/usr/bin/env python3

import argparse
import pandas as pd

JOINT_COLS = ['rt', 'mz', 'intensity', 'charge', 'quality']
SAMPLE_COLS = ['intensity']

VERSION = '1.0.0'


def main():

    args = parse_arguments()

    # ms_df = pandas.read_csv(args.input, sep=args.delim_in)
    # print(ms_df)

    cons_lines = get_consensus_lines(args.input)
    cons_fp = '{}.cons_only'.format(args.input)
    generate_consensus_intermediate(cons_fp, cons_lines)

    ms_df = pd.read_csv(cons_fp, sep=args.delim_in)
    print(ms_df.head())

    sample_numbers = get_sample_numbers(ms_df)
    normalyzer_df = setup_normalyzer_df(ms_df, sample_numbers, JOINT_COLS, SAMPLE_COLS)

    print("Writing dataframe with shape {}, {} samples to {}"
          .format(normalyzer_df.shape, len(sample_numbers), args.output))
    normalyzer_df.to_csv(args.output, sep=args.delim_out, header=False, index=False, na_rep='NA')
    print("Done!")


def get_consensus_lines(in_fp):

    cons_lines = list()
    with open(in_fp) as in_fh:
        for line in in_fh:
            line = line.rstrip()
            if line.startswith('CONSENSUS') or line.startswith('#CONSENSUS'):
                cons_lines.append(line)
    return cons_lines


def generate_consensus_intermediate(cons_fp, cons_lines):

    with open(cons_fp, 'w') as out_fh:
        for line in cons_lines:
            print(line, file=out_fh)


def get_sample_numbers(cons_df):

    all_cols = list(cons_df.columns)
    sample_nbrs = set([int(col.split('_')[-1]) for col in all_cols[1:] if not col.endswith('cf')])
    return sample_nbrs


def setup_normalyzer_df(cons_df, sample_numbers, annot_cols, sample_cols, nan_fill='NA'):

    target_sample_cols = list()
    for sample in sample_numbers:
        for col in sample_cols:
            target_sample_cols.append('{attr}_{sample}'.format(attr=col, sample=sample))

    target_annot_cols = list()
    for col in annot_cols:
        target_annot_cols.append('{}_{}'.format(col, 'cf'))

    print('TODO: Generalize sample setup to read replicate group from file')
    sample_head = [-1] + [0] * (len(target_annot_cols)-1) + [1] * len(target_sample_cols)
    sample_head_str = [str(e) for e in sample_head]

    normalyzer_vals = pd.DataFrame(cons_df[target_annot_cols + target_sample_cols].fillna(nan_fill).applymap(str))
    # normalyzer_vals_na = normalyzer_vals

    headers = pd.DataFrame([sample_head_str, list(normalyzer_vals.columns)])
    headers.columns = normalyzer_vals.columns
    normalyzer_df = pd.concat([headers, normalyzer_vals])

    return normalyzer_df


def parse_arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', help='OpenMS TSV report', required=True)
    parser.add_argument('-o', '--output', help='Normalyzer formatted report', required=True)
    parser.add_argument('--delim_in', default='\t')
    parser.add_argument('--delim_out', default='\t')
    parser.add_argument('-v', '--version', action='version', version=VERSION)

    return parser.parse_args()


if __name__ == '__main__':
    main()
