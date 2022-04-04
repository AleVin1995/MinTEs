#!/usr/bin/env python3

import numpy as np
import pandas as pd
import os
import optparse

def main():
    # create OptionParser object
    parser = optparse.OptionParser()
      
    # add options
    parser.add_option('-i', '--input',
                    dest = 'input',
                    type = 'string', 
                    help = 'path to directory containing cell line files of \
                            Bayes factors estimated by BAGEL')
    parser.add_option('-o', '--output',
                    dest = 'output',
                    type = 'string', 
                    help = 'path to save output file')
    
    (options, args) = parser.parse_args()
    
    # load input files
    cell_lines = np.char.add(options.input + '/', os.listdir(options.input))
    counter = 0
    
    for cell in cell_lines:
        if os.path.exists(cell):
            if counter == 0:
                dataset = pd.read_csv(cell, sep='\t')
                dataset = dataset.iloc[:,0:3]
            else:
                tmp = pd.read_csv(cell, sep='\t')
                tmp = tmp.iloc[:,0:3]

                dataset = pd.merge(dataset, tmp, how='outer')
        
            counter += 1
    
    if counter > 0:
        dataset.to_csv(options.output, sep='\t', index=False)


if __name__ == '__main__':
    main()