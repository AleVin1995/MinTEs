#!/usr/bin/env python3

import pandas as pd
import sys
import optparse


def main():
    # create OptionParser object
    parser = optparse.OptionParser()
      
    # add options
    parser.add_option('-i', '--input',
                    dest = 'input',
                    type = 'string', 
                    help = 'path to sgRNA-level fold-change input file')
    parser.add_option('-p', '--path',
                    dest = 'path',
                    type = 'string', 
                    help = 'path to save output files')

    (options, args) = parser.parse_args()
    
    options.input = options.input.split(' ')
    options.path = options.path.split(' ')

    for infile, outpath in zip(options.input, options.path):
        # load input file
        dataset = pd.read_csv(infile, sep='\t')
        n_cols = len(dataset.columns)

        for col_idx in range(2, n_cols):
            cell_name = dataset.columns[col_idx]
            single_cell = dataset.iloc[:,[0,1,col_idx]]
            single_cell.to_csv(outpath + '/' + cell_name + '.tsv', sep='\t', index=False)


if __name__ == '__main__':
   main()