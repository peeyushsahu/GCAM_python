__author__ = 'peeyush'
import pandas as pd
import numpy as np
import sys, os
import sompy as SOM
import SignificanceTesting as scale

def SOMclustering(Data, pheno_data, path, foldDifference, iteration = 100, gridSize=10):
    print ('Running SOM clustering')
    newDataDF = pd.DataFrame()
    pheno_groups = pheno_data.groupby('phenotype')
    #print pheno_data
    for pheno, sam in pheno_groups:
        newDataDF[pheno] = Data[list(sam['sample'].map(str.strip))].mean(axis=1)
    print 'Size of dataframe before filtering:', newDataDF.shape
    newDataDF = newDataDF[newDataDF.sum(axis=1) > 10]
    print 'Size of dataframe after filtering:', newDataDF.shape
    newDataDF['Symbol'] = newDataDF.index.str.lower()
    newDataDF.index = range(0, len(newDataDF))
    newDataDF.to_csv(os.path.join(path, 'joinedexpr.csv'), sep='\t')
    Data = newDataDF.drop(['Symbol'], axis=1)
    Data = np.array(Data) ##np.log2(Data)

    ## Grid size in rows and columns
    msz0 = gridSize
    msz1 = gridSize
    #Put this if you are updating the sompy codes otherwise simply remove it
    #reload(sys.modules['sompy'])
    sm = SOM.SOM('sm', Data, mapsize=[msz0, msz1], norm_method='var', initmethod='pca')
    sm.train(n_job=1, shared_memory='no', verbose='final', trainlen=iteration)
    sm.view_map(text_size=9, save='Yes', save_dir=os.path.join(path, '2d_plot.png'))
    #sm.hit_map(path)
    #labels = sm.cluster(method='Kmeans', n_clusters=9)
    #cents = sm.hit_map_cluster_number(path)
    cents = sm.hit_map_cluster_number(path, Data)
    #print cents[:10]
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
    #print 'Number of neurons occupied', len(df_dict)
    for ind in range(0, len(neurons)):
        df_dict[neurons[ind]] = df_dict[neurons[ind]].append(newDataDF.iloc[ind], ignore_index=True)

    return cluster_choose(df_dict, path, foldDifference=foldDifference), newDataDF


def cluster_choose(df_dict, path, foldDifference=1.9):
    clusterplot(df_dict, path)
    gene_names = []
    for df in df_dict.values():
        fpkms = []
        for a in df.mean(axis=0):
            fpkms.append(a)
        if max(fpkms)/min(fpkms) > foldDifference:
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
    '''
    Plots all the cluster of self organizing map.
    :param df_dict:
    :param path:
    :return:
    '''
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


def exprdf4plot(significanceDF, exprdata, phenodata, method, path=None, control=None, clusterSize=20):
    '''
    This methods creates a dictionary of celltypes and enriched gene expression dataframes. This dataframe is further
     analysed for the regression cofficient.
    '''
    exprdata.index = exprdata['Symbol']
    pheno = set(phenodata['phenotype'])
    expr = {}
    gene2cell_group = significanceDF.cellgenedf.groupby('CellType')
    for k, v in significanceDF.sigCelltypedf.iterrows():
        if v['P-val'] <= 0.001 and v['genecluster'] > clusterSize:
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
    if len(expr) == 0:
        raise ValueError("No significant gene for fraction comparison.")

    scaleSig = significanceDF.sigCelltypedf[(significanceDF.sigCelltypedf['P-val'] <= 0.001) & (significanceDF.sigCelltypedf['genecluster'] > clusterSize)]
    scaleSig_range = (min(scaleSig['genecluster']), max(scaleSig['genecluster']))
    for k,v in scaleSig.iterrows():
        scaleSig.loc[k, 'genecluster_scale'] = scale.scale(v['genecluster'], scaleSig_range, (1, 10))
    scaleSig.index = scaleSig['celltype']
    #print (scaleSig)
    return coffi4exprdf(expr, significanceDF, path, scaleSig, control=control, method = method)


def coffi4exprdf(expr, significanceDF, path, scaleSig, control=None, method='nuSVR'):
    '''
    Calculate cofficient and prepare dataframe for ploting.
    '''
    import statsmodels.formula.api as smf
    from sklearn.svm import NuSVR

    sigCelltypedf = significanceDF.sigCelltypedf
    sigCelltypedf = sigCelltypedf.sort(['P-val'], ascending=True)
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
    for cell in expr.keys():
        plotDataframe.loc[cell, control] = 1*scaleSig.loc[cell,'genecluster_scale']

    if method == 'one2all':
        print 'Using Multivariate regression.....'
        formula = control + ' ~ ' + ' + '.join(con)
        #print formula
        for cell in expr.keys():
            #print expr.get(cell)
            #print formula
            lm = smf.ols(formula=formula, data=expr.get(cell)).fit()
            #print lm.params
            for samp, coff in lm.params.iteritems():
                #print cell, samp, coff, scaleSig.loc[cell,'genecluster_scale']
                if coff < 0: coff = 0
                if not samp == 'Intercept':
                    plotDataframe.loc[cell, samp] = coff * scaleSig.loc[cell,'genecluster_scale']

    if method == 'one2one':
        print 'Using Univeriate regression.....'
        for cell in expr.keys():
            for condition in con:
                formula = control + ' ~ ' + condition
                #print formula
                lm = smf.ols(formula=formula, data=expr.get(cell)).fit()
                coffi = lm.params[condition]
                if lm.params[condition] < 0: coffi = 0
                plotDataframe.loc[cell, condition] = coffi * scaleSig.loc[cell,'genecluster_scale']   # plotDataframe.T

    if method == 'nuSVR':
        print 'Using nu-Support Vector Regression.....'
        for cell in expr.keys():
            data = expr.get(cell)
            target = data[control]
            trainingdf = data[con]
            #print trainingdf
            rsqr=0; nu=0.50
            for NU in [0.25, 0.50, 0.75]:
                clf = NuSVR(C=1.0, kernel='linear', nu=NU)
                clf.fit(trainingdf, target)
                rSqrd = clf.score(trainingdf,target)
                if rsqr < rSqrd:
                    rsqr = rSqrd
                    nu = NU
            print ('nu parameter:',nu)
            clf = NuSVR(C=1.0, kernel='linear', nu=nu)
            svf = clf.fit(trainingdf, target)
            coffi = svf.coef_
            #print 'Fit: ', coffi, 'rSqrd', rSqrd
            for sample in con:
                coffInd = con.index(sample)
                coff = coffi[0][coffInd]
                #print cell, sample, coff, scaleSig.loc[cell,'genecluster_scale']
                if coff < 0: coff = 0
                plotDataframe.loc[cell, sample] = coff * scaleSig.loc[cell,'genecluster_scale']
    #print plotDataframe
    if not path is None:
        plotDataframe.to_csv(os.path.join(path,'GCAM_cellexpr_sig.txt'), sep='\t')
    return plotDataframe


def plot_expressionvseignificance(path, plotdf):
    import matplotlib
    matplotlib.use('Agg')
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