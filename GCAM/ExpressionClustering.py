__author__ = 'peeyush'
import pandas as pd
import numpy as np
import sys, os
from GCAM import sompy as SOM
from GCAM import SignificanceTesting as scale
from GCAM import plots

def SOMclustering(Data, pheno_data, path, foldDifference, iteration = 100, gridSize=10, normalize=True):
    print('Running SOM clustering')
    newDataDF = pd.DataFrame()
    pheno_groups = pheno_data.groupby('phenotype')
    #print(pheno_data)
    for pheno in pheno_data['phenotype']:
        newDataDF[pheno] = Data[list(pheno_groups.get_group(pheno)['sample'].map(str.strip))].mean(axis=1)
    '''
    for pheno, sam in pheno_groups:
        newDataDF[pheno] = Data[list(sam['sample'].map(str.strip))].mean(axis=1)
    '''
    ## print 'Size of dataframe before filtering:', newDataDF.shape
    newDataDF = newDataDF[newDataDF.sum(axis=1) > 10]
    newDataDF.to_csv(os.path.join(path, 'joinedexpr.csv'), sep='\t')
    if normalize:
        ## normalizing df
        normnewDataDF = pd.DataFrame(index=newDataDF.index)
        newDataDFcolumns = newDataDF.columns
        for col in newDataDFcolumns:
            normnewDataDF[col] = (newDataDF[col] - newDataDF[col].mean())/newDataDF[col].std()
        normnewDataDF['Symbol'] = normnewDataDF.index.to_series().str.lower()
        normnewDataDF.index = range(0, len(normnewDataDF))
        normnewDataDF.to_csv(os.path.join(path, 'normjoinedexpr.txt'), sep='\t')
        Data = normnewDataDF.drop(['Symbol'], axis=1)
    else:
        normnewDataDF = newDataDF
        Data = normnewDataDF
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
    neurons = cents[:,2]
    #print 'Data points', len(neurons)
    nof_nuron_cluster = set(neurons)
    df_dict = {i:pd.DataFrame() for i in nof_nuron_cluster}
    #print 'Number of neurons occupied', len(df_dict)
    rotation = ['|', '/', '-', '\\']
    i=0
    for ind in range(0, len(neurons)):
        sys.stdout.write("\rSelecting differential genes %s" % rotation[i])
        sys.stdout.flush()
        if i == 3: i=0
        else:i+=1
        df_dict[neurons[ind]] = df_dict[neurons[ind]].append(normnewDataDF.iloc[ind], ignore_index=True)
    geneList = []
    dfGene = 2001
    while dfGene > 2000:
        print('FoldDifferencec',foldDifference)
        geneList = cluster_choose(df_dict, path, foldDifference=foldDifference)
        dfGene = len(geneList)
        with open(os.path.join(path,'parameter.txt'), 'a') as myfile:
            myfile.write("FoldChange to filter differentially expressed genes: "+str(foldDifference)+"\n")
            myfile.write("DE genes selected for analysis: "+str(dfGene)+"\n")
        foldDifference = foldDifference + 1
        myfile.close()
    return geneList, normnewDataDF


def cluster_choose(df_dict, path, foldDifference=2):
    clusterplot(df_dict, path)
    gene_names = []
    df_fold = {}
    #for df in df_dict.values():
    for key, df in df_dict.items():
        fpkms = df.mean(axis=0)
        df_fold[key] = abs(max(fpkms)/min(fpkms))
    for key in sorted(df_fold, key=df_fold.get, reverse=True):
        if df_fold[key] > foldDifference:
            for name in df_dict.get(key)['Symbol']:
                if not str(name).lower() == 'nan':
                    gene_names.append(str(name).lower())
    print('Genes are selected for analysis:', str(len(gene_names)))
    if len(gene_names) < 50:
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
    for i, v in df_dict.items():
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


def find_overlapGenefronSigcell(sigCelltypedf, cellgenedf, clustersize, args):
    '''
    This function will give common significant genes from each celltype by camparing different celltype classes.
    to one celltype.
    :param significanceDF:
    :return:
    '''
    sigCelltype = sigCelltypedf
    sigGene = cellgenedf
    cellType = []
    sigCelltype = sigCelltype.sort_values('p-val', ascending=True)
    for k, v in sigCelltype.iterrows():
        if v['genecluster'] > clustersize: # and v['p-val'] <= 0.05
            cellType.append(v['celltype'])

    sigGene_group = sigGene.groupby('celltype')
    combination = []
    genelist4expr = {}
    for cell, df in sigGene_group:
        if cell in cellType:
            #print cell, len(df)
            for cell1, df1 in sigGene_group:
                if cell1 != cell and cell1 in cellType:
                    #print cell1, len(df1)
                    if not any([cell1, cell] == x for x in combination):
                        combination.append([cell, cell1])
                        in_data1_data2 = list(set(df['gene']).intersection(df1['gene']))
                        for gene in in_data1_data2:
                            genelist4expr.setdefault(cell, []).append(gene)
                            genelist4expr.setdefault(cell1, []).append(gene)
    #print combination
    for k, v in genelist4expr.items():
        genelist4expr[k] = list(set(v))
    #print genelist4expr
    return geneList4expr(sigCelltypedf, cellgenedf, clustersize, args)
    #return find_sigCelltype4overlapGene(genelist4expr, sigCelltypedf, cellgenedf, clustersize, args)


def geneList4expr(sigCelltypedf, cellgenedf, clustersize, args):
    sigCelltype = sigCelltypedf
    sigGene = cellgenedf
    cellType = []
    sigCelltype = sigCelltype.sort_values('p-val', ascending=True)
    for k, v in sigCelltype.iterrows():
        if v['genecluster'] > clustersize: # and v['p-val'] <= 0.05
            cellType.append(v['celltype'])
    genelist4expr = dict((k, []) for k in cellType)
    sigGene_group = sigGene.groupby('gene')
    for gene, df in sigGene_group:
        df = df.sort_values('p-val', ascending=True)
        for k, v in df.iterrows():
            if v['celltype'] in cellType:
                genelist4expr[v['celltype']].append(v['gene'])
                break
    print(args.outdir)
    myfile = open(os.path.join(args.outdir, 'celltypeSign_new.txt'), 'w')
    myfile.write('celltype'+'\t'+'genes')
    for k, v in genelist4expr.items():
        myfile.write('\n'+k+'\t'+','.join(v))
    myfile.close()
    return genelist4expr




def find_sigCelltype4overlapGene(genelist4expr, sigCelltypedf, cellgenedf, clusterSize, args):
    '''
    Find for which celltype a overlap gene is significant
    :return:
    '''
    sigCelltype = sigCelltypedf
    sigGene = cellgenedf
    cellType = genelist4expr.keys()
    dict4sortedGenes = {}
    gene2cell_group = sigGene.groupby('celltype')
    # print genelist4expr
    # Check in celltype sig. table for sig cell type.
    if len(sigCelltype[sigCelltype['genecluster'] > clusterSize]) < 3:
        sys.exit('Warning: Less than 3 cell-types are above cluster size. '
                 'Please choose a lower threshold for Cluster size for significant cell-type.')
    for k, v in sigCelltype.iterrows():
        if v['genecluster'] > clusterSize:  # and v['p-val'] <= 0.05
            # print ('Celltype above cluster size:', v['celltype'], v['p-val'], v['genecluster'])
            df = gene2cell_group.get_group(v['celltype'])
            # Take the sig cell type and get genes associated to it.
            for key, val in df.iterrows():
                # check if the associated genes are in the overlaping list
                if val['gene'] not in genelist4expr[v['celltype']]:
                    dict4sortedGenes.setdefault(v['celltype'], []).append(val['gene'])
                else:
                    # gene in a overapping list test which celltype has the most significant pval for the gene.
                    celltype = None; enrichment = 0
                    for sigcell in cellType:
                        gene_ind = gene2cell_group.get_group(sigcell)
                        if val['gene'] in list(gene_ind['gene']):
                            ind = gene_ind['gene'][gene_ind['gene'] == val['gene']].index[0]
                            if gene_ind.loc[ind, 'enrichment'] > enrichment:
                                enrichment = gene_ind.loc[ind, 'enrichment']
                                celltype = sigcell
                                # print sigcell, pval
                    if celltype is not None: dict4sortedGenes.setdefault(celltype, []).append(val['gene'])

    for cell, sigGenes in dict4sortedGenes.iteritems():
        dict4sortedGenes[cell] = list(set(sigGenes))
        # print cell, len(sigGenes)
    print("Genes/celltype after filtering:")
    #for key, cell in dict4sortedGenes.iteritems():
        #print(key, len(cell))
    myfile = open(os.path.join(args.outdir, 'celltypeSign_old.txt'), 'w')
    myfile.write('celltype'+'\t'+'genes')
    for k, v in dict4sortedGenes.iteritems():
        myfile.write('\n'+k+'\t'+','.join(v))
    myfile.close()
    return dict4sortedGenes


def exprdf4plot(significanceDF, exprdata, phenodata, args, path=None, control=None, clusterSize=20):
    '''
    This methods creates a dictionary of celltypes and enriched gene expression dataframes. This dataframe is further
     analysed for the regression cofficient.
    '''
    from collections import OrderedDict

    exprdata.index = exprdata['Symbol']
    #pheno = set(phenodata['phenotype'])
    pheno = phenodata['phenotype']
    userCelltype = args.selectCelltypes
    expr = {}
    sigCelltypedf = significanceDF.sigCelltypedf
    cellgenedf = significanceDF.cellgenedf
    if userCelltype is not None:
        userCelltype = [x.strip() for x in userCelltype.split(',')]
        sigCelltypedf = sigCelltypedf[sigCelltypedf['celltype'].isin(userCelltype)]
        cellgenedf = cellgenedf[cellgenedf['celltype'].isin(userCelltype)]
    #print sigCelltypedf
    '''
    if not args.key_celltype_list:
        sigCelltypedf = sigCelltypedf[sigCelltypedf['genecluster'] > clusterSize]
    '''
    if userCelltype is None:
        sigCelltypedf = significanceDF.sigCelltypedf
        sigCelltypedf = sigCelltypedf[sigCelltypedf['genecluster'] > clusterSize]

    # if args.remOverlapping:
    print('Removing overlapping genes')
    gene2cellCluster = find_overlapGenefronSigcell(sigCelltypedf, cellgenedf, clusterSize, args)
    cellexpr_dict = OrderedDict()
    for celltype, genelist in gene2cellCluster.items():
        if len(genelist) > 10:  # Minimum genes per celltype for comparision
            #print celltype
            for i in pheno:
                geneexpr_dict = {}
                for gene in genelist:
                    geneexpr_dict[gene] = exprdata.loc[gene, i]
                cellexpr_dict[i] = geneexpr_dict
            expr[celltype] = pd.DataFrame(cellexpr_dict)
            #expr[celltype].to_csv(os.path.join(path, celltype+'_GCAM_cellexpr.txt'), sep='\t')
            #print expr
    if len(expr) == 0:
        raise ValueError("No significant gene for fraction comparison.")
    if args.controlsample is None:
        print("mean as control.")
        return mean_coffi4exprdf(expr, path, args)
    return coffi4exprdf(expr, path, args, control=control)


def coffi4exprdf(expr, path, args, control=None):
    '''
    Calculate cofficient and prepare dataframe for ploting.
    '''
    import statsmodels.formula.api as smf
    from sklearn.svm import NuSVR
    con = []
    plotDataframe = pd.DataFrame(index=expr.keys())
    for sample in expr.keys():
        for cont in expr.get(sample).columns:
            if not cont == control:
                con.append(cont)
        break
    #print(control, con)
    print('Using nuSVR.....')
    for cell in expr.keys():
        data = expr.get(cell)
        target = data[control]
        #print con
        for col in con:
            trainingdf = data[[col]]
            rsqr=0; nu=0.50
            for NU in [0.25, 0.50, 0.75]:
                clf = NuSVR(C=1.0, kernel='linear', nu=NU)
                clf.fit(trainingdf, target)
                rSqrd = clf.score(trainingdf,target)
                if rsqr < rSqrd:
                    rsqr = rSqrd
                    nu = NU
            #print ('nu parameter:',nu)
            clf = NuSVR(C=1.0, kernel='linear', nu=nu)
            svf = clf.fit(trainingdf, target)
            coffi = svf.coef_
            if coffi < 0: coffi = 0
            plotDataframe.loc[cell, col+'_vs_'+control] = coffi #* scaleSig.loc[cell,'genecluster_scale']
    #print plotDataframe
    if not path is None:
        plotDataframe.to_csv(os.path.join(path, 'GCAM_cellexpr_sig.txt'), sep='\t')
        for column in plotDataframe.columns:
            plotDataframe.loc[:,column] = (plotDataframe.loc[:,column]/sum(plotDataframe.loc[:,column]))
    ## plotting stacked bar plot
    plots.heatmap_Sigcelltype(plotDataframe.T, path)
    plots.stack_barplot(args, plotDataframe, path, name=control, method='exprbased')
    #return plotDataframe

def mean_coffi4exprdf(expr, path, args):
    '''
    Calculate cofficient and prepare dataframe for ploting.
    '''
    from sklearn.svm import NuSVR
    con = []
    plotDataframe = pd.DataFrame(index=expr.keys())
    for sample in expr.keys():
        for cont in expr.get(sample).columns:
            con.append(cont)
        break
    print('nuSVR.....')
    for cell in expr.keys():
        #print(cell)
        data = expr.get(cell)
        target = data.sum(axis=1)
        #data.to_csv(os.path.join(path, 'GCAM_reg_exprdata.txt'), sep='\t')
        for col in con:
            trainingdf = data[[col]]
            rsqr=0; nu=0.50
            for NU in [0.25, 0.50, 0.75]:
                clf = NuSVR(C=1.0, kernel='linear', nu=NU)
                clf.fit(trainingdf, target)
                rSqrd = clf.score(trainingdf,target)
                #print("Value of nu:", NU," rSqrd: ", rSqrd)
                if rsqr < rSqrd:
                    rsqr = rSqrd
                    nu = NU
            #print ('nu parameter:',nu)
            clf = NuSVR(C=1.0, kernel='linear', nu=nu)
            svf = clf.fit(trainingdf, target)
            coffi = svf.coef_
            if coffi < 0 or rsqr < 0: coffi = 0
            plotDataframe.loc[cell, col] = coffi #* scaleSig.loc[cell,'genecluster_scale']
        df2exprsig = pd.concat([target, data], axis=1)
        #df2exprsig.to_csv(os.path.join(path, cell+'GCAM_cellexpr_sig.txt'), sep='\t')

    #print plotDataframe
    if not path is None:
        plotDataframe.to_csv(os.path.join(path, 'GCAM_cellexpr_sig.txt'), sep='\t')
        for column in plotDataframe.columns:
            plotDataframe.loc[:,column] = (plotDataframe.loc[:,column]/sum(plotDataframe.loc[:,column]))
    ## plotting stacked bar plot
    plots.heatmap_Sigcelltype(plotDataframe.T, path)
    plots.stack_barplot(args, plotDataframe, path, name='_', method='exprbased')
    #return plotDataframe


def plot_expressionvseignificance(path, plotdf):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    import math
    plotdf = plotdf[plotdf['p-val'] < 0.05]
    print('plotting significance plot')
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
