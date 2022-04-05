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
                    help = 'path to Bayes factor dataset')
    parser.add_option('-e', '--essential', 
                    dest = 'essential',
                    type = 'string', 
                    help = 'path to reference essential genes file')
    parser.add_option('-n', '--nonessential', 
                    dest = 'nonessential',
                    type = 'string', 
                    help = 'path to reference nonessential genes file')
    parser.add_option('--fdr', 
                    dest = 'false_discovery',
                    type = 'float', 
                    help = 'correction at X% fdr rate')
    parser.add_option('-o', '--output',
                    dest = 'output',
                    type = 'string', 
                    help = 'path to save output file')
    
    (options, args) = parser.parse_args()

    # load input file
    dataset = pd.read_csv(options.input, sep='\t')

    essential = pd.read_csv(options.essential, sep='\t')
    essential = essential.to_numpy().flatten()

    nonessential = pd.read_csv(options.nonessential, sep='\t')
    nonessential = nonessential.to_numpy().flatten()

    for cell_idx in range(2, len(dataset.columns)):
        thr = 1-options.false_discovery
        cell_line = dataset.columns[cell_idx]

        bf_ess = dataset.loc[dataset.iloc[:,1].isin(essential), cell_line]
        bf_ess = bf_ess[~bf_ess.isna()]
        bf_ess = np.array(bf_ess)

        bf_non = dataset.loc[dataset.iloc[:,1].isin(nonessential), cell_line]
        bf_non = bf_non[~bf_non.isna()]
        bf_non = np.array(bf_non)

        label = np.concatenate(([1] * len(bf_ess), [0] * len(bf_non)))
        bf_train = np.concatenate((bf_ess, bf_non))

        bf_ind_sort = np.argsort(bf_train)
        bf_train_sort = bf_train[bf_ind_sort]

        labelcum = np.cumsum(np.flip(label[bf_ind_sort]))
        ppv = labelcum / np.arange(1, len(labelcum)+1)
        ppv = np.flip(ppv)

        thr_pos = np.where(ppv >= thr)[0]

        # In case not possible to scale at X% FDR
        # It may happen for low-quality screens
        if len(thr_pos) == 0:
            thr_pos = np.where(ppv == max(ppv))[0]

        thr_bf = np.min(bf_train_sort[thr_pos])
        dataset.iloc[:, cell_idx] -= thr_bf
    
    dataset.to_csv(options.output, sep='\t', index=False)

if __name__ == '__main__':
    main()