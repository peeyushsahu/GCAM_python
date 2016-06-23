from GCAM.plots import HiearchicalHeatmap
__author__ = 'peeyush'
import pandas as pd
import scipy.stats as stats
import os
import logging, timeit
from GCAM import FilesFolders as read
from GCAM import plots



class SignificanceObject():
    '''
    This Object will hold dataframes used and created in analysis
    '''
    def __init__(self, occurrencedf, binom_prob, resourcepath, outdir, heatmapdf=None):
        self.outdir = outdir
        self.resourcepath = resourcepath
        self.occurrencedf = occurrencedf
        self.binom_prob = binom_prob
        self.heatmapdf = heatmapdf
        self.filheatmapdf = None
        self.pvaldf = None
        self.adjpvaldf = None
        self.cellgenedf = None
        self.sigCelltypedf = None
        self.binom_pval_df = None

    def heatmapdf_create(self, thres):
        '''
        This function will generate df for HeatMap by scaling the occurrencedf
        '''
        occurrencedf = self.occurrencedf
        transposedf = pd.DataFrame.transpose(occurrencedf)
        #print transposedf.head()
        scaled_df = scale_dataframe(transposedf)
        #print scaled_df.head()
        #print occurrencedf.index
        #print scaled_df
        scaled_df.columns = occurrencedf.index
        scaled_df = scaled_df.set_index(occurrencedf.columns)
        self.heatmapdf = scaled_df
        self.filter_heatmapdf(thres=thres)

    def filter_heatmapdf(self, thres=(20,20)):
        '''
        This method filters rows and columns with sum < 1 in HeatMapdf
        '''
        df = self.heatmapdf
        self.filheatmapdf = df.loc[df.sum(1) > thres[0], df.sum(0) > thres[1]]

    def plot_heatmap(self, path):
        '''
        This method will plot the HeatMap dataframe. Using package HclustHeatMap.
        :param HiearchicalHeatmap:
        :param df:
        :param path:
        :return:
        '''
        hclustHeatmap = HiearchicalHeatmap()
        hclustHeatmap.frame = self.filheatmapdf
        hclustHeatmap.path = os.path.sep.join([path, 'GCAM_heatMap.svg'])
        fig, axm, axcb, cb = hclustHeatmap.plot()

    def filter_occuDf(self):
        '''
        Filter occurrence df for removing gene wid less than 5 celltype tags
        :return:
        '''
        occuDf = self.occurrencedf
        Columns=[]
        for k, v in occuDf.iteritems():
            #print k, v.sum()
            if v.sum() < 5:
                Columns.append(k)
        #print(len(Columns))
        self.occurrencedf = occuDf.drop(Columns, axis=1)
        #print(self.occurrencedf.shape)


    def fisher_occurrence_test(self):
        '''
        This method will calculate significance of celltypes per gene using their occurrence.
        Statistical test used is Fisher Exact Test
        '''
        fsstart = timeit.default_timer()
        self.filter_occuDf()
        occu_df = self.occurrencedf
        binom_prob = self.binom_prob
        pvaldf = pd.DataFrame()
        binom_pvaldf = pd.DataFrame()
        adjpvaldf = pd.DataFrame()
        enrichmentdf = pd.DataFrame()
        matsum = occu_df.sum().sum()
        #print occu_df.head()
        for k, v in occu_df.iterrows():
            key = v.keys()
            rowsum = v.sum()
            for i in range(0, v.shape[0]):
                value = v[i]
                if not value == 0:
                    colsum = occu_df[[i]].sum()[0] - value
                    rsum = rowsum - value
                    #print rowsum, value, colsum
                    enrichment = float(value)/occu_df[[i]].sum()[0]
                    if value != 0:
                        ## Fisher p-value is calcualted but not put in the table
                        oddsratio, pval = stats.fisher_exact([[value, colsum], [rsum, matsum-(value+rsum+colsum)]],
                                                             alternative='greater')
                        celltype = k
                        index = binom_prob[binom_prob['celltype'] == celltype].index.tolist()
                        b_pval = stats.binom.sf(value, colsum+value, binom_prob.iloc[index[0], 3])

                    else:
                        pval = 1
                        b_pval = 1
                    pvaldf.loc[k, key[i]] = pval
                    binom_pvaldf.loc[k, key[i]] = b_pval
                    enrichmentdf.loc[k, key[i]] = enrichment
                else:
                    binom_pvaldf.loc[k, key[i]] = 1
                    pvaldf.loc[k, key[i]] = 1
                    enrichmentdf.loc[k, key[i]] = 0
        ## adj p-val calcualtion
        #print(binom_pvaldf.head())
        for k, v in binom_pvaldf.iterrows():
            for i in range(0, v.shape[0]):
                key = v.keys()
                value = v[i]
                if value != 1:
                    #print(key[i])
                    #print(len(binom_pvaldf[binom_pvaldf[key[i]] < 0.05]))
                    sigCelltype = len(binom_pvaldf[binom_pvaldf[key[i]] < 0.05])
                    #print(sigCelltype)
                    if value < 0.05/sigCelltype:
                        adjpvaldf.loc[k, key[i]] = value*sigCelltype
                    else:
                        adjpvaldf.loc[k, key[i]] = 1
                else:
                    adjpvaldf.loc[k, key[i]] = 1
        #print binom_pvaldf
        self.binom_pval_df = binom_pvaldf
        self.pvaldf = pvaldf
        self.adjpvaldf = adjpvaldf
        fsstop = timeit.default_timer()
        logging.info('TC in sig occur test:'+str(fsstop-fsstart)+'sec')
        self.celltype_overrepresntation_list(enrichmentdf)  #### def()

    def celltype_overrepresntation_list(self, enrichmentdf):
        '''
        This method will save the result of significance in one DF.
        '''
        significance = 1
        column = ['celltype', 'gene', 'enrichment', 'p-val', 'FDR']
        cellgenedf = pd.DataFrame()
        for celltype, v in self.binom_pval_df.iterrows():
            for gene, pval in v.iteritems():
                if pval < significance:
                    cellgenedf = cellgenedf.append(pd.Series([celltype, gene, enrichmentdf.loc[celltype, gene], pval,
                                                         self.adjpvaldf.loc[celltype, gene]]), ignore_index=True)
        #print cellgenedf.head(10)
        cellgenedf.columns = column
        print('cellgenedf shape:', cellgenedf.shape)
        #cellgenedf = self.filter_df(cellgenedf)
        self.cellgenedf = cellgenedf
        print('cellgenedf shape after:', cellgenedf.shape)
        #self.filter_cellgenedf()  # Filter single cell multigene enrihment
        self.fisher_significant_celltypes()

    def fisher_significant_celltypes(self):
        '''
        This method will test the combined significance of celltype in the data and help predicts
        its association with user given data.
        '''
        fsstart = timeit.default_timer()
        sigcelltype = pd.DataFrame()
        #print self.cellgenedf
        cellgroup = self.cellgenedf.groupby(self.cellgenedf['celltype'])
        cellgenedf = self.cellgenedf
        c = len(cellgenedf[cellgenedf['p-val'] <= 0.05])
        d = len(cellgenedf) #[cellgenedf['P-val'] < 0.5]
        #print(c,d)
        for celltype, val in cellgroup:
            #print celltype
            if len(val[val['p-val'] <= 0.001]) > 1:
                #print val
                a = len(val[val['p-val'] <= 0.001])
                b = len(val) - a
                cc = c - a
                dd = d - (a+b+cc)
                #print a, ':', b, ':', cc, ':', dd, c, d
                oddsRatio, p = stats.fisher_exact([[a, b], [cc, dd]])
                #print 'celltype:'+celltype, a, p
                sigcelltype = sigcelltype.append(pd.Series([celltype, a, p]), ignore_index=True)
        sigcelltype.columns = ['celltype', 'genecluster', 'p-val']
        length = len(sigcelltype[sigcelltype['p-val'] <= 0.05])
        # Significant celltype check
        if length < 1:
            raise ValueError('No siginificant celltypes.')

        for k, v in sigcelltype.iterrows():
            if v['p-val'] < 0.05/length:
                sigcelltype.loc[k, 'FDR'] = v['p-val']*length
            else:
                sigcelltype.loc[k, 'FDR'] = 1
        self.sigCelltypedf = sigcelltype
        fsstop = timeit.default_timer()
        logging.info('TC in sig celltype test:'+str(fsstop-fsstart)+'sec')


    def data4radarplot(self):
        '''
        Preparing data for radar plot
        '''
        import numpy as np
        tissue = []
        genes = []
        sigCelltypedf = self.sigCelltypedf
        cell2tissue = read.cell2tissue_DB(self.resourcepath)
        #print(cell2tissue.head(3))
        cell2tissue['genes'] = 0
        for k, v in sigCelltypedf.iterrows():
            if v['celltype'] in list(cell2tissue.index):
                cell2tissue.loc[v['celltype'], 'genes'] = v['genecluster']

        norm_cell2gene = cell2tissue['tissue'].value_counts()
        #print(norm_cell2gene)
        df_cell_gr = cell2tissue.groupby('tissue')
        for cell in df_cell_gr.groups:
            #print(cell)
            df = df_cell_gr.get_group(cell)
            tissue.append(cell)
            genes.append((sum(df['genes'])*1.)/norm_cell2gene.loc[cell])
        gene2plot = np.divide(genes, sum(genes))
        #print(gene2plot, tissue, self.outdir)
        plots.plot_radar(gene2plot, tissue, self.outdir)



def scale(val, src, dst):
    '''
    This returns scaled value
    :param val: value to be scaled
    :param src: min and max of values to be scaled
    :param dst: range to scale
    :return:
    '''
    return (((val - src[0]) * (dst[1] - dst[0])) / (src[1] - src[0])) + dst[0]


def scale_dataframe(df):
    '''
    This will scale a Pandas.DataFrame with every row min-max
    :param df: DataFrame
    :return: Scaled DataFrame
    '''
    new_max = 100
    new_min = 0
    list_of_rows = []
    #print 'Process: scaling of dataframe'
    for r, v in df.iterrows():
        #print v.sum()
        rows = []
        for val in v:
            #print val
            old_min = min(v)
            old_max = max(v)
            if not isinstance(val, str):
                #print val
                rows.append(scale(val, (old_min, old_max), (new_min, new_max)))
                #rows.append(float(val)/v.sum())
        list_of_rows.append(rows)
    return pd.DataFrame(data=list_of_rows)
