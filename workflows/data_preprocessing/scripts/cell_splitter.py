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
    
    #options.input = options.input.split(' ')
    #options.path = options.path.split(' ')

    # load input file
    dataset = pd.read_csv(options.input, sep='\t')
    n_cols = len(dataset.columns)
    
    for col_idx in range(2, n_cols):
        cell_idx = col_idx-1
        single_cell = dataset.iloc[:,[0,1,col_idx]]
        single_cell.to_csv(options.path + '/cell_line_' + str(cell_idx) + '.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()