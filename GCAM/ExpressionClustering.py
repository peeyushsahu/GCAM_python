__author__ = 'peeyush'
import pandas as pd
import time as time
import numpy as np
from matplotlib import pyplot as plt
#pd.__version__
import sys, os
import sompy as SOM

def SOMclustering(Data, pheno_data, path):
    #Data = pd.read_csv('/home/peeyush/Downloads/GSE72502_PBMC_RPKMs_1.csv', sep=',', header=0, index_col=0)
    #pheno_data = pd.read_csv('/home/peeyush/Downloads/phenotype.csv', sep=',', header=0)
    #path = '/home/peeyush/Downloads'
    newDataDF = pd.DataFrame()
    pheno_groups = pheno_data.groupby('phenotype')
    print pheno_data
    for pheno, sam in pheno_groups:
        newDataDF[pheno] = Data[list(sam['sample'].map(str.strip))].mean(axis=1)
    newDataDF['Symbol'] = newDataDF.index.str.lower()
    newDataDF.index = range(0, len(newDataDF))
    newDataDF.to_csv(os.path.join(path, 'joinedexpr.csv'), sep='\t')
    Data = newDataDF.drop(['Symbol'], axis=1)

    ## The critical factor which increases the computational time, but mostly the memory problem is the size of SOM (i.e. msz0,msz1),
    ## other wise the training data will be parallelized
    ## This is your selected map size
    msz0 = 10
    msz1 = 10
    #This is a random data set, but in general it is assumed that you have your own data set as a numpy ndarray
    #Data = np.random.rand(10*1000,20)
    Data = np.array(Data)

    #Put this if you are updating the sompy codes otherwise simply remove it
    #reload(sys.modules['sompy'])
    sm = SOM.SOM('sm', Data, mapsize=[msz0, msz1], norm_method='var', initmethod='pca')
    sm.train(n_job=1, shared_memory='no', verbose='final', trainlen=1000)
    sm.view_map(text_size=9, save='Yes', save_dir=os.path.join(path, '2d_plot.png'))
    sm.hit_map(path)
    labels = sm.cluster(method='Kmeans', n_clusters=9)
    #cents = sm.hit_map_cluster_number(path)
    cents = sm.hit_map_cluster_number(path, Data)
    print cents[:10]
    '''
    clustering = labels[cents[:,2]]
    nof_cluster = set(clustering)
    print ('Number of clusters:',len(nof_cluster))
    df_dict = {i:pd.DataFrame() for i in nof_cluster}
    for ind in range(0, len(clustering)):
        df_dict[clustering[ind]] = df_dict[clustering[ind]].append(newDataDF.iloc[ind], ignore_index=True)
    '''
    neurons = cents[:,2]
    #print 'Data points', len(neurons)
    nof_nuron_cluster = set(neurons)
    df_dict = {i:pd.DataFrame() for i in nof_nuron_cluster}
    print 'Number of neurons occupied', len(df_dict)
    for ind in range(0, len(neurons)):
        df_dict[neurons[ind]] = df_dict[neurons[ind]].append(newDataDF.iloc[ind], ignore_index=True)

    return cluster_choose(df_dict, path), newDataDF


def cluster_choose(df_dict, path):
    clusterplot(df_dict, path)
    gene_names = []
    for df in df_dict.values():
        fpkms = []
        for a in df.mean(axis=0):
            fpkms.append(a)
        if max(fpkms)/min(fpkms) > 6:
            #print max(fpkms), min(fpkms)
            for name in df['Symbol']:
                if not str(name).lower() == 'nan':
                    #print name
                    gene_names.append(str(name).lower())
    print 'Genes are selected for analysis:', len(gene_names)
    if len(gene_names) < 20:
        raise ValueError("Not enough differential genes.")
    return gene_names


def clusterplot(df_dict, path):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(nrows=10, ncols=10)
    fig.set_size_inches(18.5, 10.5)
    row = 0
    col = 0
    for i, v in df_dict.iteritems():
        v.boxplot(ax=axes[row,col], fontsize=7)
        axes[row,col].set_title('Size:'+str(len(v)), fontsize=10)
        row += 1
        #print row, col
        if row == 10:
            row = 0
            col += 1
        if i == 99: break
    fig.tight_layout()
    fig.savefig(os.path.join(path, 'GCAM_cluster_expression.png'))
    fig.clf()


def exprdf4plot(significanceDF, exprdata, phenodata, path = None, control=None):
    '''
    This methods creates a dictionary of celltypes and enriched gene expression dataframes. This dataframe is further
     analysed for the regression cofficient.
    '''
    exprdata.index = exprdata['Symbol']
    pheno = set(phenodata['phenotype'])
    expr = {}
    gene2cell_group = significanceDF.cellgenedf.groupby('CellType')
    for k, v in significanceDF.sigCelltypedf.iterrows():
        if v['P-val'] <= 0.001:
            cellexpr_dict = {}
            df = gene2cell_group.get_group(v['celltype'])
            for i in pheno:
                geneexpr_dict = {}
                for key, val in df.iterrows():
                    if val['P-val'] <= 0.05:
                        geneexpr_dict[val['Genes']] = exprdata.loc[val['Genes'], i]
                cellexpr_dict[i] = geneexpr_dict
            expr[v['celltype']] = pd.DataFrame(cellexpr_dict)
    #print expr
    return coffi4exprdf(expr, significanceDF, path, control)


def coffi4exprdf(expr, significanceDF, path, control=None):
    import statsmodels.formula.api as smf
    sigCelltypedf = significanceDF.sigCelltypedf
    sigCelltypedf.index = sigCelltypedf['celltype']
    plotDataframe = pd.DataFrame(index=expr.keys())
    expr.keys()
    con = []
    for sample in expr.keys():
        #print sample
        #print expr.get(sample).columns
        for cont in expr.get(sample).columns:
            if not cont == control:
                con.append(cont)
        break
    '''
    formula = control + ' ~ ' + '+'.join(con)
    #print formula
    for cell in expr.keys():
        #print expr.get(cell)
        lm = smf.ols(formula=formula, data=expr.get(cell)).fit()
        print lm.params
        for samp, coff in lm.params.iteritems():
            if coff < 0: coff = 0
            if not samp == 'Intercept': plotDataframe.loc[cell, samp] = coff
    '''
    for cell in expr.keys():
        for condition in con:
            formula = control + ' ~ ' + condition
            lm = smf.ols(formula=formula, data=expr.get(cell)).fit()
            plotDataframe.loc[cell, condition] = lm.params[condition]
    # plotDataframe.T
    if not path is None:
        plotDataframe.to_csv(os.path.join(path,'GCAM_cellexpr_sig.txt'), sep='\t')
    return plotDataframe


def plot_expressionvseignificance(path, plotdf):
    import matplotlib.pyplot as plt
    import numpy as np
    import math
    plotdf = plotdf[plotdf['p-val'] < 0.05]
    print ('plotting significance plot')
    plotdf = plotdf.sort(['P-val'], ascending=True)
    l = plotdf['genecluster'].tolist()
    t = plotdf['P-val'].tolist()
    s = plotdf['relative_expression'].tolist()
    name = plotdf['celltype'].tolist()
    area = [(math.log(x, 10) * -15) for x in t]
    color = np.random.random(len(t))

    plt.scatter(l, s, s=area, c=color, alpha=0.5)
    plt.grid(True, linestyle=':', color='black')
    # draw a thick red hline at y=0 that spans the xrange
    h = plt.axhline(linewidth=1, color='r', linestyle='--')

    # draw a default vline at x=1 that spans the yrange
    h = plt.axvline(linewidth=1, color='r', x=0.001, linestyle='--')
    for i in range(0, len(t)):
        plt.annotate(name[i], xy=(l[i], s[i]), xycoords='data',
            xytext=(-10,1), textcoords='offset points',
            ha='center', va='bottom',
            bbox=dict(boxstyle='round, pad=0.2', fc='yellow', alpha=0.2),
            fontsize=8)
## Plot legend
    l1 = plt.scatter([],[], s=50, c='gray', alpha=0.5)
    l2 = plt.scatter([],[], s=200, c='gray', alpha=0.5)
    labels = ["less significant", "highly significant"]
    plt.legend([l1, l2], labels, ncol=2, frameon=True, fontsize=8,
    handlelength=2, loc = 4, borderpad = 0.5,
    handletextpad=1, scatterpoints=1)

    plt.tick_params(axis='both', labelsize=8)
    plt.xlim(0, max(l)+20)
    plt.ylim(min(s)-2,max(s)+20)
    plt.title('Relative Cell-type expression plot', fontsize=14)
    plt.xlabel('Gene cluster size', fontsize=12)
    plt.ylabel('Average Fold Change', fontsize=12)
    #plt.tight_layout()
    plt.savefig(path + os.path.sep + 'GCAM_celltype_VS_expresiion.png')
    plt.clf()