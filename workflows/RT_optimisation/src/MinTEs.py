#!/usr/bin/env python3

from os import supports_bytes_environ
import numpy as np
import pandas as pd
import math
import scipy.stats as stats
from sklearn.linear_model import LinearRegression
import optparse


def load_data(fold_change, bayes_factor, essential_genes, non_essential_genes):
    """
    Load data

    Args:
        fold_change: dataframe (tsv format) containing depletion fold-change scores.
                     First column must contain single-guide RNA (sgRNA) IDs.
                     Second column must contain the gene symbols targeted by the 
                     corresponding sgRNA.
                     All the remaining columns represent CRISPR-Cas9 screens (i.e., 
                     cell line IDs), containing fold-change scores derived from the
                     sgRNA i targeting gene j in screen c.
        bayes_factor: same format as per the fold_change argument.  
                      First column must contain single-guide RNA (sgRNA) IDs.
                      Second column must contain the gene symbols targeted by the 
                      corresponding sgRNA.
                      All the remaining columns represent CRISPR-Cas9 screens (i.e., 
                      cell line IDs), containing bayes factors derived from the
                      corresponding screen in the fold_change dataframe using BAGEL2
                      with user-provided sets of full reference essential and non-
                      essential genes.
        essential_genes: dataframe (tsv format) containing the full reference set of
                         essential genes. 
                         First column must contain gene symbol (other columns are not
                         necessary).
        nonessential_genes: dataframe (tsv format) containing the full reference set 
                            of essential genes.
                            First column must contain gene symbol (other columns are 
                            not necessary).
    """
    
    fc = pd.read_table(fold_change)
    bf = pd.read_table(bayes_factor)

    coreEss = pd.read_table(essential_genes)
    coreEss = coreEss.iloc[np.in1d(coreEss.iloc[:,0],fc.iloc[:,1]),0]
    coreEss = coreEss.to_numpy().flatten()

    nonEss = pd.read_table(non_essential_genes)
    nonEss = nonEss.iloc[np.in1d(nonEss.iloc[:,0],fc.iloc[:,1]),0]
    nonEss = nonEss.to_numpy().flatten()

    return(fc, bf, coreEss, nonEss)


def rounding(tot, sub):
    """
    Rounding BAGEL-outputted bayes factor array, obtained with reduce templates, to 
    the number of decimals present in the total array (default 4)

    Args:
        tot: array of bayes factors using the full reference gene templates.
        sub: array of bayes factors using the reduced gene templates at a specific
             according to the parameters (CRISPR-Cas9 library, reference gene sets
             and percentage of subsampling).
    """

    tot = str(tot)
    tot = tot[1:(len(tot) - 1)].split(' ')

    round_dec = max([len(i.split('.')[1]) for i in tot if i != ''])
    sub = np.round(sub, round_dec)
    return sub


def testing(fold_change, bayes_factor, essential_genes, non_essential_genes,
                sub_ess, sub_non, scaling):
    fc, bf, coreEss, nonEss = fold_change, bayes_factor, essential_genes, non_essential_genes

    """
    Testing performances of pre-defined reduce templates vs full gene templates

    Args:
        fold_change: dataframe (tsv format) containing depletion fold-change scores.
                     First column must contain single-guide RNA (sgRNA) IDs.
                     Second column must contain the gene symbols targeted by the 
                     corresponding sgRNA.
                     All the remaining columns represent CRISPR-Cas9 screens (i.e., 
                     cell line IDs), containing fold-change scores derived from the
                     sgRNA i targeting gene j in screen c.
        bayes_factor: same format as per the fold_change argument.  
                      First column must contain single-guide RNA (sgRNA) IDs.
                      Second column must contain the gene symbols targeted by the 
                      corresponding sgRNA.
                      All the remaining columns represent CRISPR-Cas9 screens (i.e., 
                      cell line IDs), containing bayes factors derived from the
                      corresponding screen in the fold_change dataframe using BAGEL2
                      with user-provided sets of full reference essential and non-
                      essential genes.
        essential_genes: dataframe (tsv format) containing the full reference set of
                         essential genes. 
                         First column must contain gene symbol (other columns are not
                         necessary).
        nonessential_genes: dataframe (tsv format) containing the full reference set 
                            of essential genes.
                            First column must contain gene symbol (other columns are 
                            not necessary).
        sub_ess: dataframe (tsv format) containing the reduced reference set of 
                 essential genes. 
        sub_non: dataframe (tsv format) containing the reduced reference set of 
                 nonessential genes. 
        scaling: boolean parameter, whether to scale computed bayes factor to 5%
                 false discovery rate (default is True).
    """

    # 1st column must contain single-guide RNA (sgRNA) IDs
    # 2nd column must contain the gene symbols targeted by the corresponding sgRNA
    # All the remaining columns represent CRISPR-Cas9 screens (i.e., cell line IDs)
    fc = fc.sort_values(fc.columns[0])
    bf = bf.sort_values(bf.columns[0])
    
    # consider only common cell lines
    cells = fc.columns[2:].to_numpy()
    cells = np.intersect1d(cells, bf.columns)
    
    # testing phase
    JScell = np.zeros(len(cells))

    train_ess = fc[fc.iloc[:, 1].isin(sub_ess)]
    train_non = fc[fc.iloc[:, 1].isin(sub_non)]
    
    for i in range(0, len(cells)):
        cellLine = cells[i]
        
        train_guides = fc[cellLine]
        train_guides = train_guides[~train_guides.isna()]
        train_guides = np.array(train_guides)

        train_ess_fc = train_ess[cellLine]
        train_ess_fc = train_ess_fc[~train_ess_fc.isna()]
        train_ess_fc = np.array(train_ess_fc)

        train_non_fc = train_non[cellLine]
        train_non_fc = train_non_fc[~train_non_fc.isna()]
        train_non_fc = np.array(train_non_fc)

        bf_tot = bf[cellLine]
        bf_tot = bf_tot[~bf_tot.isna()]
        bf_tot = bf_tot.to_numpy().flatten()

        bfcell = calculate_bayes_factors(train_guides, train_ess_fc, train_non_fc, len(sub_non), scaling)

        # compute Jaccard similarity for each cell line
        if type(bfcell) is np.ndarray:
            bfcell = rounding(bf_tot, bfcell)
            JScell[i] = len(np.where((bf_tot > 0.0) & (bfcell > 0.0))[0]) / len(np.where((bf_tot > 0.0) | (bfcell > 0.0))[0])
        else:
            JScell[i] = 0.0
    
    JS_test = np.mean(JScell)
    JS_test = np.round(JS_test, 4)

    # prepare output dictionary with all summary data
    top = max(len(coreEss), len(nonEss), len(cells), len(JScell), len(sub_ess), len(sub_non))
    Xnan = [float("nan")]

    res = dict()
    res["JS_test"] = [JS_test] + Xnan * (top - 1)
    res["scaling"] = [scaling] + Xnan * (top - 1)
    res["cells"] = cells.tolist() + Xnan * (top - len(cells))
    res["coreEss"] = coreEss.tolist() + Xnan * (top - len(coreEss))
    res["nonEss"] = nonEss.tolist() + Xnan * (top - len(nonEss))
    res["combo_ess"] = sub_ess.tolist() + Xnan * (top - len(sub_ess))
    res["combo_non"] = sub_non.tolist() + Xnan * (top - len(sub_non))
    return(res)


def greedy(fold_change, bayes_factor, essential_genes, non_essential_genes,
                subsample, step, fraction_train, seed, scaling, sub_rate, const):
    fc, bf, coreEss, nonEss = fold_change, bayes_factor, essential_genes, non_essential_genes
    frTrain = fraction_train

    """
    Greedy-search strategy implemented by MinTEs to select optimised reduced templates according
    to user-defined parameters (i.e., CRISPR-Cas9 library, reference gene sets and percentage of 
    subsampling)

    Args:
        fold_change: dataframe (tsv format) containing depletion fold-change scores.
                     First column must contain single-guide RNA (sgRNA) IDs.
                     Second column must contain the gene symbols targeted by the 
                     corresponding sgRNA.
                     All the remaining columns represent CRISPR-Cas9 screens (i.e., 
                     cell line IDs), containing fold-change scores derived from the
                     sgRNA i targeting gene j in screen c.
        bayes_factor: same format as per the fold_change argument.  
                      First column must contain single-guide RNA (sgRNA) IDs.
                      Second column must contain the gene symbols targeted by the 
                      corresponding sgRNA.
                      All the remaining columns represent CRISPR-Cas9 screens (i.e., 
                      cell line IDs), containing bayes factors derived from the
                      corresponding screen in the fold_change dataframe using BAGEL2
                      with user-provided sets of full reference essential and non-
                      essential genes.
        essential_genes: dataframe (tsv format) containing the full reference set of
                         essential genes. 
                         First column must contain gene symbol (other columns are not
                         necessary).
        nonessential_genes: dataframe (tsv format) containing the full reference set 
                            of essential genes.
                            First column must contain gene symbol (other columns are 
                            not necessary).
        subsample: percentage of subsampling (default is 0.05) in order to obtain
                   reduced gene templates.
        step: nº iterations without improvement (default is 10) before decreasing the 
              substitution rate (in case of a dynamic substitution rate).
        fraction_train: Fraction of cell lines used for training (0,1] in case of 
                        multiple screens (i.e., more than one cell line screened).
        seed: method to initialise random number generator (optional).
        scaling: boolean parameter, whether to scale computed bayes factor to 5%
                 false discovery rate (default is True).
        sub_rate: substitution rate applied between the subset of reference genes and
                  the hold-out genes (default is 1).
        const: boolean parameter, whether to apply a static substitution rate (default
               is False).
    """

    np.random.seed(seed)

    # assume 1st column contains sgRNA IDs
    # 2nd column the corresponding gene IDs
    # the rest the cell line IDs
    fc = fc.sort_values(fc.columns[0])
    bf = bf.sort_values(bf.columns[0])
    
    # consider only common cell lines
    cells = fc.columns[2:].to_numpy()
    common_cols = np.intersect1d(cells, bf.columns)
    
    count = 0
    keep = True
    
    # training and test cell lines
    if frTrain == 1.0 or len(common_cols) == 1:
        cells_train = common_cols
        cells_test = common_cols
    else:
        np.random.shuffle(common_cols)
        cells_train = common_cols[:round((len(common_cols) * frTrain))]
        cells_test = np.setdiff1d(common_cols, cells_train)

    frEss = len(coreEss) / (len(coreEss) + len(nonEss))
    takein = round((len(coreEss) + len(nonEss)) * subsample)
    Ein = min(round(takein * frEss), len(coreEss))
    Nin = takein - Ein

    JS = np.zeros(takein * max(step, 1))
    sr_rec = np.zeros(takein * max(step, 1))
    
    # training phase
    while keep:
        JScell = np.zeros(len(cells_train))

        if count == 0:
            combo_ess = np.random.choice(coreEss, size = Ein, replace = False)
            combo_non = np.random.choice(nonEss, size = Nin, replace = False)

            sub_rate_ess = int(Ein * sub_rate)
            sub_rate_non = int(Ein * sub_rate)
        else:
            combo_ess_old = combo_ess.copy()
            combo_non_old = combo_non.copy()

            sub_rate_ess = max(int(Ein * sub_rate), 1)
            sub_rate_non = max(int(Nin * sub_rate), 1)

            np.random.shuffle(combo_ess)
            np.random.shuffle(combo_non)

            ess_out = np.setdiff1d(coreEss, combo_ess[sub_rate_ess:])
            non_out = np.setdiff1d(nonEss, combo_non[sub_rate_non:])

            np.random.shuffle(ess_out)
            np.random.shuffle(non_out)

            combo_ess[:sub_rate_ess] = ess_out[:sub_rate_ess]
            combo_non[:sub_rate_non] = non_out[:sub_rate_non]

        train_ess = fc[fc.iloc[:, 1].isin(combo_ess)]
        train_non = fc[fc.iloc[:, 1].isin(combo_non)]
        
        for i in range(0, len(cells_train)):
            cellLine = cells_train[i]
            
            train_guides = fc[cellLine]
            train_guides = train_guides[~train_guides.isna()]
            train_guides = np.array(train_guides)

            train_ess_fc = train_ess[cellLine]
            train_ess_fc = train_ess_fc[~train_ess_fc.isna()]
            train_ess_fc = np.array(train_ess_fc)

            train_non_fc = train_non[cellLine]
            train_non_fc = train_non_fc[~train_non_fc.isna()]
            train_non_fc = np.array(train_non_fc)

            bf_tot = bf[cellLine]
            bf_tot = bf_tot[~bf_tot.isna()]
            bf_tot = bf_tot.to_numpy().flatten()

            bfcell = calculate_bayes_factors(train_guides, train_ess_fc, train_non_fc, len(combo_non), scaling)

            # compute Jaccard similarity for each cell line
            if type(bfcell) is np.ndarray:
                bfcell = rounding(bf_tot, bfcell)
                JScell[i] = len(np.where((bf_tot > 0.0) & (bfcell > 0.0))[0]) / len(np.where((bf_tot > 0.0) | (bfcell > 0.0))[0])
            else:
                JScell[i] = 0.0

        JS[count] = np.mean(JScell)
        JS[count] = np.round(JS[count], 4)
        sr_rec[count] = np.round(sub_rate, 4)

        # extend JS and sr_rec arrays if needed
        if JS[-1] != 0:
            JS = np.concatenate((JS, np.zeros(takein * step)))
            sr_rec = np.concatenate((sr_rec, np.zeros(takein * step)))
        
        # keep the old combo if the JS is worse
        if count > 0:
            if JS[count] <= JS[count - 1]:
                JS[count] = JS[count - 1]
                combo_ess, combo_non = combo_ess_old, combo_non_old
        
        # stop algorithm if there is no improvement for count consecutive steps
        if count >= step:
            if JS[count] == JS[count - step] and const == True:
                keep = False
            if JS[count] == JS[count - step] and sr_rec[count] == sr_rec[count - step]:
                if sub_rate_ess > 1 or sub_rate_non > 1:
                    sub_rate *= np.exp(-1)
                else:
                    keep = False

        count += 1

    JS = JS[JS != 0]
    JS_train = JS[-1]

    sr_rec = sr_rec[sr_rec != 0]
    
    # test phase
    if frTrain == 1.0 or len(cells) == 1:
        JS_test = JS_train
    else:
        JScell = np.zeros(len(cells_test))
        test_ess = fc[fc.iloc[:,1].isin(combo_ess)]
        test_non = fc[fc.iloc[:,1].isin(combo_non)]

        for i in range(0, len(cells_test)):
            cellLine = cells_test[i]

            test_guides = fc[cellLine]
            test_guides = test_guides[~test_guides.isna()]
            test_guides = np.array(test_guides)

            test_ess_fc = test_ess[cellLine]
            test_ess_fc = test_ess_fc[~test_ess_fc.isna()]
            test_ess_fc = np.array(test_ess_fc)

            test_non_fc = test_non[cellLine]
            test_non_fc = test_non_fc[~test_non_fc.isna()]
            test_non_fc = np.array(test_non_fc)

            bf_tot = bf[cellLine]
            bf_tot = bf_tot[~bf_tot.isna()]
            bf_tot = bf_tot.to_numpy().flatten()

            bfcell = calculate_bayes_factors(test_guides, test_ess_fc, test_non_fc, len(combo_non), scaling)

            # compute Jaccard similarity for each cell line
            if type(bfcell) is np.ndarray:
                JScell[i] = len(np.where((bf_tot > 0.0) & (bfcell > 0.0))[0]) / len(np.where((bf_tot > 0.0) | (bfcell > 0.0))[0])
            else:
                JScell[i] = 0.0

        JS_test = np.mean(JScell)
        JS_test = np.round(JS_test, 4)
    
    # create output dictionary with all summary data
    top = max(len(coreEss), len(nonEss), len(cells_train), len(cells_test), len(JS))
    Xnan = [float("nan")]

    res = dict()
    res["JS_train"] = [JS_train] + Xnan * (top - 1)
    res["JS_test"] = [JS_test] + Xnan * (top - 1)
    res["History"] = JS.tolist() + Xnan * (top - len(JS))
    res["step"] = [step] + Xnan * (top - 1)
    res["seed"] = [seed] + Xnan * (top - 1)
    res["sub_rate"] = sr_rec.tolist() + Xnan * (top - len(sr_rec))
    res["scaling"] = [scaling] + Xnan * (top - 1)
    res["fraction_train"] = [frTrain] + Xnan * (top - 1)
    res["subsample"] = [subsample] + Xnan * (top - 1)
    res["cells_train"] = cells_train.tolist() + Xnan * (top - len(cells_train))
    res["cells_test"] = cells_test.tolist() + Xnan * (top - len(cells_test))
    res["coreEss"] = coreEss.tolist() + Xnan * (top - len(coreEss))
    res["nonEss"] = nonEss.tolist() + Xnan * (top - len(nonEss))
    res["combo_ess"] = combo_ess.tolist() + Xnan * (top - len(combo_ess))
    res["combo_non"] = combo_non.tolist() + Xnan * (top - len(combo_non))
    return(res)


def calculate_bayes_factors(train_guides, train_ess, train_non, n_nonEss, scaling = False):
    """
    Compute bayes factors from depletion fold-change score in a CRISPR-Cas9 screen

    Args:
        train_guides: depletion fold-change scores of guides targeting reference genes (either full
                      reduced gene sets)
        train_ess: character vector of guides targeting essential reference genes (either full or
                   reduced gene sets).
        train_non: character vector of guides targeting nonessential reference genes (either full or
                   reduced gene sets).
        n_nonEss: nº of nonessential genes.
        scaling: boolean parameter, whether to scale computed bayes factor to 5% false discovery 
                 rate (default is False in this function).
    """

    FC_THRESH = 2 ** (-1.1535 * np.log(n_nonEss + 13.324) + 0.7728)
    
    # gaussian kernel fit onto essential and nonessential training discrete distributions
    kess = stats.gaussian_kde(train_ess)
    knon = stats.gaussian_kde(train_non)
    
    # define empirical upper and lower bounds within which to calculate BF = f(fold change)
    x = np.arange(-10, 2, 0.01)
    nonfitx = knon.evaluate(x)

    # define lower bound empirical fold change threshold:  minimum FC np.where knon is above threshold
    f = np.where(nonfitx > FC_THRESH)
    xmin = np.around(min(x[f]*100))/100.0

    # define upper bound empirical fold change threshold:  minimum value of log2(ess/non)
    subx = np.arange(xmin, max(x[f]), 0.01)
    logratio_sample = np.log2(kess.evaluate(subx) / knon.evaluate(subx))
    f = np.where(logratio_sample == logratio_sample.min())
    xmax = np.around(subx[f]*100)/100.0
    xmax = xmax[0]

    # round foldchanges to nearest 0.01
    # precalculate logratios and build lookup table (for speed)
    span = np.arange(xmin, xmax + 0.01, 0.01)
    logratio_lookup = np.log2(kess.evaluate(span) / knon.evaluate(span))

    # linear interpolation
    foldchange = train_guides[np.where((train_guides >= xmin) & (train_guides <= xmax))]

    if len(foldchange) == 0:
        return None

    foldchange = np.around(foldchange * 100)
    testx = foldchange / 100.0
    coord = foldchange - np.min(foldchange)
    testy = logratio_lookup[coord.astype(int)]

    slope, intercept, r_value, p_value, std_err = stats.linregress(np.array(testx), np.array(testy))
    
    # Bayes factor calculation
    bf = slope * train_guides + intercept

    # Bayes factor scaling
    if scaling:
        bf_ess = slope * train_ess + intercept
        bf_non = slope * train_non + intercept

        label = np.concatenate(([1] * len(bf_ess), [0] * len(bf_non)))
        bf_train = np.concatenate((bf_ess, bf_non))

        bf_ind_sort = np.argsort(bf_train)
        bf_train_sort = bf_train[bf_ind_sort]

        labelcum = np.cumsum(np.flip(label[bf_ind_sort]))
        ppv = labelcum / np.arange(1, len(labelcum)+1)
        ppv = np.flip(ppv)

        thr_pos = np.where(ppv >= 0.95)[0]

        # In case not possible to scale at 5% FDR
        # It may happen for low-quality screens
        if len(thr_pos) == 0:
            thr_pos = np.where(ppv == max(ppv))[0]

        thr_bf = np.min(bf_train_sort[thr_pos])
        bf -= thr_bf

    return(bf)


def main():
    # create OptionParser object
    parser = optparse.OptionParser()
      
    # add options
    parser.add_option("-e", "--essential", 
                    dest = "essential",
                    type = "string", 
                    help = "path to reference essential genes file")
    parser.add_option("-n", "--nonessential", 
                    dest = "nonessential",
                    type = "string", 
                    help = "path to reference nonessential genes file")
    parser.add_option("-f", "--fold_change", 
                    dest = "fold_change",
                    type = "string", 
                    help = "specify path to fold change file")
    parser.add_option("-b", "--bayes_factor",
                    dest = "bayes_factor",
                    type = "string", 
                    help = "specify path to reference bayesian factor file")
    parser.add_option("-o", "--output",
                    dest = "output",
                    type = "string", 
                    help = "specify output file name")
    parser.add_option("-s", "--seed",
                    dest = "seed",
                    type = "int",
                    default = np.random.randint(1e9), 
                    help = "specify seed (Optional)")
    parser.add_option("--subsample",
                    dest = "subsample",
                    type = "float",
                    default = 0.05, 
                    help = "Percentage of subsampling (0,1]")
    parser.add_option("--fr_train",
                    dest = "fraction_train",
                    type = "float",
                    default = 1.0, 
                    help = "Fraction of cell lines used for training (0,1]")
    parser.add_option("--step",
                    dest = "step",
                    type = "int",
                    default = 10, 
                    help = "nº steps without improvement")
    parser.add_option("--scaling",
                    dest = "scaling",
                    action = "store_true",
                    default = False, 
                    help = "scaling of bayes factors")
    parser.add_option("--sub_rate",
                    dest = "sub_rate",
                    type = "float",
                    default = 1.0, 
                    help = "Initial substitution rate (0,1]")
    parser.add_option("-c", "--constant",
                    dest = "const",
                    action = "store_true",
                    default = False, 
                    help = "constant substitution rate")
    parser.add_option("--subess",
                    dest = "subess",
                    default = None, 
                    help = "pre-defined set of reduced reference essential genes")
    parser.add_option("--subnon",
                    dest = "subnon",
                    default = None, 
                    help = "pre-defined set of reduced reference nonessential genes")

    (options, args) = parser.parse_args()
    fc, bf, coreEss, nonEss = load_data(options.fold_change, options.bayes_factor,
                                            options.essential, options.nonessential)

    # In case reduce templates are provided, MinTEs will compare the performances
    # of these reduced templates in predicting fitness genes against those of the
    # full reference templates
    if options.subess != None and options.subnon != None:
        sub_ess = pd.read_table(options.subess)
        sub_ess = sub_ess.loc[sub_ess.iloc[:,0].isin(coreEss),:]
        sub_ess = sub_ess.iloc[:,0]
        sub_ess = sub_ess.to_numpy().flatten()

        sub_non = pd.read_table(options.subnon)
        sub_non = sub_non.loc[sub_non.iloc[:,0].isin(nonEss),:]
        sub_non = sub_non.iloc[:,0]
        sub_non = sub_non.to_numpy().flatten()
        res = testing(fc, bf, coreEss, nonEss, sub_ess, sub_non, 
                        options.scaling)
    
    # Execute MinTEs optimisation to search for the best combination of reduced
    # templates given user-specified parameters (i.e., CRISPR-Cas9 library, 
    # reference gene sets and percentage of subsampling)
    else:
        res = greedy(fc, bf, coreEss, nonEss, options.subsample, 
                        options.step, options.fraction_train, 
                        options.seed, options.scaling,
                        options.sub_rate, options.const)

    res = pd.DataFrame(res)
    res.to_csv(options.output, sep = '\t', index = False)


if __name__ == '__main__':
    main()
