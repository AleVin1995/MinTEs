#!/usr/bin/env python3

import pandas as pd
import optparse


def main():
    # create OptionParser object
    parser = optparse.OptionParser()
      
    # add options
    parser.add_option('-i', '--input',
                    dest = 'input',
                    type = 'string', 
                    help = 'path to sgRNA-level fold-change input file')
    parser.add_option('-o', '--output',
                    dest = 'output',
                    type = 'string',
                    help = 'path to save output files')

    (options, args) = parser.parse_args()
    
    # load input file
    dataset = pd.read_csv(options.input, sep='\t')
    max_cells = len(dataset.columns)

    for cell_idx in range(2, max_cells):
        single_cell = dataset.iloc[:,[0,1,cell_idx]]
        single_cell = single_cell.dropna()
        single_cell.to_csv(options.output + '/cell_line_' + str(cell_idx-1) + '.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()