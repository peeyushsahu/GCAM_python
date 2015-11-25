__author__ = 'peeyush'
import pandas as pd
import numpy as np
import sys, os
import sompy as SOM
import SignificanceTesting as scale
import plots

def SOMclustering(Data, pheno_data, path, foldDifference, iteration = 100, gridSize=10):
    print('Running SOM clustering')
    newDataDF = pd.DataFrame()
    pheno_groups = pheno_data.groupby('phenotype')
    #print(pheno_data)
    for pheno, sam in pheno_groups:
        newDataDF[pheno] = Data[list(sam['sample'].map(str.strip))].mean(axis=1)
    ## print 'Size of dataframe before filtering:', newDataDF.shape
    newDataDF = newDataDF[newDataDF.sum(axis=1) > 10]
    newDataDF.to_csv(os.path.join(path, 'joinedexpr.csv'), sep='\t')
    ## normalizing df
    normnewDataDF = pd.DataFrame(index=newDataDF.index)
    newDataDFcolumns = newDataDF.columns
    for col in newDataDFcolumns:
        normnewDataDF[col] = (newDataDF[col] - newDataDF[col].mean())/newDataDF[col].std()
    normnewDataDF['Symbol'] = normnewDataDF.index.str.lower()
    normnewDataDF.index = range(0, len(normnewDataDF))
    normnewDataDF.to_csv(os.path.join(path, 'normjoinedexpr.csv'), sep='\t')
    Data = normnewDataDF.drop(['Symbol'], axis=1)
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
    for ind in range(0, len(neurons)):
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
    for key, df in df_dict.iteritems():
        fpkms = df.mean(axis=0)
        df_fold[key] = abs(max(fpkms)/min(fpkms))
    for key in sorted(df_fold, key=df_fold.get, reverse=True):
        if df_fold[key] > foldDifference:
            for name in df_dict.get(key)['Symbol']:
                if not str(name).lower() == 'nan':
                    gene_names.append(str(name).lower())
    print('Genes are selected for analysis:', len(gene_names))
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


def find_overlapGenefronSigcell(sigCelltypedf, cellgenedf, clustersize):
    '''
    This function will give common significant genes from each celltype by camparing diggerent celltype classes.
    to one celltype.
    :param significanceDF:
    :return:
    '''
    sigCelltype = sigCelltypedf
    sigGene = cellgenedf
    cellType = []
    sigCelltype = sigCelltype.sort('P-val', ascending=True)
    for k, v in sigCelltype.iterrows():
        if v['P-val'] < 0.001 and v['genecluster'] > clustersize:
            cellType.append(v['celltype'])

    sigGene_group = sigGene.groupby('CellType')
    combination = []
    genelist4expr = {}
    #sigGene_group.get_group('alveolar macrophage').head(2)
    for cell, df in sigGene_group:
        if cell in cellType:
            #print cell, len(df)
            for cell1, df1 in sigGene_group:
                if cell1 != cell and cell1 in cellType:
                    #print cell1, len(df1)
                    if not any([cell1, cell] == x for x in combination):
                        combination.append([cell, cell1])
                        in_data1_data2 = list(set(df['Genes']).intersection(df1['Genes']))
                        for gene in in_data1_data2:
                            genelist4expr.setdefault(cell, []).append(gene)
                            genelist4expr.setdefault(cell1, []).append(gene)
    #print combination
    for k, v in genelist4expr.iteritems():
        #print k, len(v)
        #print v
        genelist4expr[k] = list(set(v))
    #print genelist4expr
    return find_sigCelltype4overlapGene(genelist4expr, sigCelltypedf, cellgenedf, clusterSize=clustersize)


def find_sigCelltype4overlapGene(genelist4expr, sigCelltypedf, cellgenedf, clusterSize=20):
    '''
    Find for which celltype a overlap gene is significant
    :return:
    '''
    sigCelltype = sigCelltypedf
    sigGene = cellgenedf
    cellType = genelist4expr.keys()
    dict4sortedGenes = {}
    gene2cell_group = sigGene.groupby('CellType')
    ## Check in celltype sig. table for sig cell type.
    for k, v in sigCelltype.iterrows():
        if v['P-val'] <= 0.001 and v['genecluster'] > clusterSize:
            df = gene2cell_group.get_group(v['celltype'])
            ## Take the sig cell type and get genes associated to it.
            for key, val in df.iterrows():
                if val['P-val'] <= 0.05:
                    ## check if the associated genes are in the overlaping list
                    if val['Genes'] not in genelist4expr[v['celltype']]:
                        dict4sortedGenes.setdefault(v['celltype'], []).append(val['Genes'])
                    else:
                        ## gene in a overapping list test which celltype has the most significant pval for the gene.
                        celltype = None; pval = 1
                        for sigcell in cellType:
                            #print sigcell
                            gene_ind = gene2cell_group.get_group(sigcell)
                            #print sigcell
                            if val['Genes'] in list(gene_ind['Genes']):
                                #print 'True'
                                ind = gene_ind['Genes'][gene_ind['Genes'] == val['Genes']].index[0]
                                if gene_ind.loc[ind, 'P-val'] < pval:
                                    if sigcell not in ['t lymphocyte', 'neuron', 'stem cell', 'b lymphocyte','fibroblast', 'muscle cell']: # onyl assign sig cell type if it not among these.
                                        pval = gene_ind.loc[ind, 'P-val']
                                        celltype = sigcell
                                        #print sigcell, pval
                        if celltype is not None: dict4sortedGenes.setdefault(celltype, []).append(val['Genes'])
                        #print 'Final: ',celltype, pval
    for cell, sigGenes in dict4sortedGenes.iteritems():
        dict4sortedGenes[cell] = list(set(sigGenes))
        #print cell, len(sigGenes)
    return dict4sortedGenes


def exprdf4plot(significanceDF, exprdata, phenodata, args, path=None, control=None, clusterSize=20):
    '''
    This methods creates a dictionary of celltypes and enriched gene expression dataframes. This dataframe is further
     analysed for the regression cofficient.
    '''
    exprdata.index = exprdata['Symbol']
    pheno = set(phenodata['phenotype'])
    expr = {}
    if not args.key_celltype_list:
        sigCelltypedf = significanceDF.sigCelltypedf
        sigCelltypedf = sigCelltypedf.sort('P-val', ascending=True)
        sigCelltypedf = sigCelltypedf[:15]
        #print(sigCelltypedf['celltype'])
    else:
        sigCelltypedf = significanceDF.sigCelltypedf
    cellgenedf = significanceDF.cellgenedf
    if args.remOverlapping:
        print('Removing overlapping genes')
        gene2cellCluster = find_overlapGenefronSigcell(sigCelltypedf, cellgenedf, clusterSize)
        cellexpr_dict = {}
        for celltype, genelist in gene2cellCluster.iteritems():
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

    else:
        gene2cell_group = cellgenedf.groupby('CellType')
        for k, v in sigCelltypedf.iterrows():
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

    scaleSig = sigCelltypedf[(sigCelltypedf['P-val'] <= 0.001) & (sigCelltypedf['genecluster'] > clusterSize)]
    scaleSig_range = (min(scaleSig['genecluster']), max(scaleSig['genecluster']))
    for k,v in scaleSig.iterrows():
        scaleSig.loc[k, 'genecluster_scale'] = scale.scale(v['genecluster'], scaleSig_range, (1, 10))
    scaleSig.index = scaleSig['celltype']
    #print (scaleSig)
    if args.meanAsControl:
        print("mean as control.")
        return mean_coffi4exprdf(expr, path, scaleSig, args, method = args.regMethod)
    return coffi4exprdf(expr, path, scaleSig, args, control=control, method = args.regMethod)


def coffi4exprdf(expr, path, scaleSig, args, control=None, method='nuSVR'):
    '''
    Calculate cofficient and prepare dataframe for ploting.
    '''
    import statsmodels.formula.api as smf
    from sklearn.svm import NuSVR
    '''
    samples = []
    for sample in expr.keys():
        for cont in expr.get(sample).columns:
            samples.append(cont)
        break
    '''
    #for control in samples:
    con = []
    plotDataframe = pd.DataFrame(index=expr.keys())
    for sample in expr.keys():
        for cont in expr.get(sample).columns:
            if not cont == control:
                con.append(cont)
        break
    #print(control, con)
    for cell in expr.keys():
        plotDataframe.loc[cell, 'Reference_%'] = 1#*scaleSig.loc[cell,'genecluster_scale']
    if method == 'one2all':
        print ('Using Multivariate regression.....')
        formula = control + ' ~ ' + ' + '.join(con)
        #print formula
        for cell in expr.keys():
            #print expr.get(cell)
            #print formula
            lm = smf.ols(formula=formula, data=expr.get(cell)).fit()
            #print lm.params
            for samp, coff in lm.params.iteritems():
                #print cell, samp, coff, scaleSig.loc[cell,'genecluster_scale']
                #if coff < 0: coff = 0
                if not samp == 'Intercept':
                    plotDataframe.loc[cell, samp+'_vs_reference'] = coff * scaleSig.loc[cell,'genecluster_scale']

    if method == 'one2one':
        print('Using Univeriate regression.....')
        for cell in expr.keys():
            for condition in con:
                formula = control + ' ~ ' + condition
                #print formula
                lm = smf.ols(formula=formula, data=expr.get(cell)).fit()
                coffi = lm.params[condition]
                #if lm.params[condition] < 0: coffi = 0
                plotDataframe.loc[cell, condition+'_vs_reference'] = coffi * scaleSig.loc[cell,'genecluster_scale']   # plotDataframe.T

    if method == 'nuSVR':
        print('Using nu-Support Vector Regression.....')
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
            #print ('nu parameter:',nu)
            clf = NuSVR(C=1.0, kernel='linear', nu=nu)
            svf = clf.fit(trainingdf, target)
            coffi = svf.coef_
            #print 'Fit: ', coffi, 'rSqrd', rSqrd
            #print con
            for sample in con:
                coffInd = con.index(sample)
                coff = coffi[0][coffInd]
                if coff < 0: coff = 0
                plotDataframe.loc[cell, sample+'_vs_'+control] = coff #* scaleSig.loc[cell,'genecluster_scale']

    if method == 'nuSVR-one2one':
        print('Using Univeriate nuSVR.....')
        for cell in expr.keys():
            data = expr.get(cell)
            target = data[control]
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
        for column in plotDataframe.columns:
            plotDataframe.loc[:,column] = (plotDataframe.loc[:,column]/sum(plotDataframe.loc[:,column]))
            #plotDataframe.loc[:,column] = (plotDataframe.loc[:,column]+abs(plotDataframe.loc[:,column].mean()))
        #plotDataframe.to_csv(os.path.join(path,control+'GCAM_cellexpr_sig.txt'), sep='\t')
    ## plotting stacked bar plot
    plots.heatmap_Sigcelltype(plotDataframe.T, path)
    plots.stack_barplot(plotDataframe, path, name=control, key_celltypes=args.key_celltype_list, method='exprbased')
    #return plotDataframe

def mean_coffi4exprdf(expr, path, scaleSig, args, method='nuSVR'):
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
    '''
    for cell in expr.keys():
        plotDataframe.loc[cell, 'Reference_%'] = 1#*scaleSig.loc[cell,'genecluster_scale']
    '''
    if method == 'nuSVR':
        print('Using nu-Support Vector Regression.....')
        for cell in expr.keys():
            data = expr.get(cell)
            target = data.sum(axis=1)
            data.to_csv(os.path.join(path, 'GCAM_reg_exprdata.txt'), sep='\t')
            trainingdf = data[con]
            rsqr=0; nu=0.50
            for NU in [0.25, 0.50, 0.75]:
                clf = NuSVR(C=1.0, kernel='linear', nu=NU)
                clf.fit(trainingdf, target)
                rSqrd = clf.score(trainingdf,target)
                if rsqr < rSqrd:
                    rsqr = rSqrd
                    nu = NU
            clf = NuSVR(C=1.0, kernel='linear', nu=nu)
            svf = clf.fit(trainingdf, target)
            coffi = svf.coef_
            for sample in con:
                coffInd = con.index(sample)
                coff = coffi[0][coffInd]
                if coff < 0: coff = 0
                plotDataframe.loc[cell, sample] = coff # scaleSig.loc[cell,'genecluster_scale']

    if method == 'nuSVR-one2one':
        print('Using Univeriate nuSVR.....')
        for cell in expr.keys():
            #print(cell)
            data = expr.get(cell)
            target = data.sum(axis=1)
            data.to_csv(os.path.join(path, 'GCAM_reg_exprdata.txt'), sep='\t')
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

    #print plotDataframe
    if not path is None:
        for column in plotDataframe.columns:
            plotDataframe.loc[:,column] = (plotDataframe.loc[:,column]/sum(plotDataframe.loc[:,column]))
            #plotDataframe.loc[:,column] = (plotDataframe.loc[:,column]+abs(plotDataframe.loc[:,column].mean()))
        #plotDataframe.to_csv(os.path.join(path, 'GCAM_cellexpr_sig.txt'), sep='\t')
    ## plotting stacked bar plot
    plots.heatmap_Sigcelltype(plotDataframe.T, path)
    plots.stack_barplot(plotDataframe, path, name='_', key_celltypes=args.key_celltype_list, method='exprbased')
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