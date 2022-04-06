#!/usr/bin/env python3

import numpy as np
import pandas as pd
import os
import optparse

def split_callback(option, opt, value, parser):
    setattr(parser.values, option.dest, value.split(' '))

def main():
    # create OptionParser object
    parser = optparse.OptionParser()
      
    # add options
    parser.add_option('--RT',
                    dest = 'RT',
                    type = 'string',
                    action='callback',
                    callback=split_callback,
                    help = 'path to directory containing RT sets')
    parser.add_option('--cross',
                    dest = 'cross',
                    type = 'string',
                    action='callback',
                    callback=split_callback,
                    help = 'path to directory containing cross-testings')
    parser.add_option('-o', '--output',
                    dest = 'output',
                    type = 'string', 
                    help = 'path to save output file')
    
    (options, args) = parser.parse_args()
    
    # Build ranking dataframe
    counter = 0

    for RT_folder, cross_folder in zip(options.RT, options.cross):
        percentages = os.listdir(RT_folder)

        for perc in percentages:
            RT_files = os.listdir(os.path.join(RT_folder, perc))
            cross_files = os.listdir(os.path.join(cross_folder, perc))

            for RT_test, cross_test in zip(RT_files, cross_files):
                RT_path = os.path.join(RT_folder, perc, RT_test)

                RT_summary = pd.read_csv(RT_path, sep='\t')
                cross_summary = pd.read_csv(os.path.join(cross_folder, perc, cross_test), sep='\t')

                items = RT_path.split('/')
                resources_idx = np.where(np.array(items) == 'resources')[0][0]
                items = items[resources_idx:]

                essential_gene_set = items[2]
                project = items[3].split('_RT')[0]
                subsampling = items[4]

                iteration = items[5].split('_')[1]
                iteration = int(iteration.split('.')[0])

                JS_train = float(RT_summary['JS_train'].dropna()[0])
                JS_test = float(RT_summary['JS_test'].dropna()[0])
                JS_cross = float(cross_summary['JS_test'].dropna()[0])

                if counter == 0:
                    n_rows = len(options.RT) * len(percentages) * len(RT_files)
                    idx = [float("nan") for i in range(1, n_rows+1)]
                    column_names = ['Essential_gene_set', 'Project', 'Subsampling',
                                    'Iteration', 'JS_train', 'JS_test', 'JS_cross', 
                                    'lib_dep_RT_score', 'lib_indep_RT_score']

                    ranking_df = pd.DataFrame(dict((column, idx) for column in column_names))
                
                ranking_df.loc[counter, 'Essential_gene_set'] = essential_gene_set
                ranking_df.loc[counter, 'Project'] = project
                ranking_df.loc[counter, 'Subsampling'] = subsampling
                ranking_df.loc[counter, 'Iteration'] = iteration
                ranking_df.loc[counter, 'JS_train'] = JS_train
                ranking_df.loc[counter, 'JS_test'] = JS_test
                ranking_df.loc[counter, 'JS_cross'] = JS_cross
                ranking_df.loc[counter, 'lib_dep_RT_score'] = JS_train*0.2 + JS_test*0.8
                ranking_df.loc[counter, 'lib_indep_RT_score'] = JS_train*0.1 + JS_test*0.4 + JS_cross*0.5

                counter+=1

    ranking_df.to_csv(options.output, sep='\t', index=False)


if __name__ == '__main__':
    main()