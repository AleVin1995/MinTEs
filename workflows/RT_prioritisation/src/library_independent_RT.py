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
                    help = 'path to file ranking RT sets')
    parser.add_option('-o', '--output',
                    dest = 'output',
                    type = 'string', 
                    help = 'path to directory to save output files')
    
    (options, args) = parser.parse_args()

    # load input file
    ranking_df = pd.read_csv(options.input, sep='\t')

    # selection top library specific RT
    idx = ranking_df.groupby(['Essential_gene_set', 'Subsampling'])['lib_indep_RT_score'].transform(max) == ranking_df['lib_indep_RT_score']
    ranking_df = ranking_df[idx]
    ranking_df = ranking_df.reset_index()

    for set_idx in range(0, len(ranking_df)):
        essential_gene_set = ranking_df.loc[set_idx, 'Essential_gene_set']
        project = ranking_df.loc[set_idx, 'Project']
        subsampling = str(ranking_df.loc[set_idx, 'Subsampling'])
        iteration = str(int(ranking_df.loc[set_idx, 'Iteration']))
        
        RT_ess = pd.DataFrame()
        RT_noness = pd.DataFrame()
        
        test_file = pd.read_csv(os.path.join('resources/BF', essential_gene_set, project + '_RT', subsampling, 'test_' + iteration + '.txt'), sep='\t')
        
        RT_ess['GENE'] = test_file['combo_ess'].dropna()
        RT_noness['GENE'] = test_file['combo_non'].dropna()

        RT_ess.to_csv(os.path.join(options.output, essential_gene_set + '_' + subsampling + '_ess.tsv'), sep='\t', index=False)
        RT_noness.to_csv(os.path.join(options.output, essential_gene_set + '_' + subsampling + '_noness.tsv'), sep='\t', index=False)


if __name__ == '__main__':
    main()