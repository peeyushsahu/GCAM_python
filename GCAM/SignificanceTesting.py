from GCAM.plots import HiearchicalHeatmap
__author__ = 'peeyush'
import pandas as pd
import scipy.stats as stats
import os


class SignificanceObject():
    '''
    This Object will hold dataframes used and created in analysis
    '''
    def __init__(self, occurrencedf, heatmapdf=None):
        self.occurrencedf = occurrencedf
        self.heatmapdf = heatmapdf
        self.filheatmapdf = None
        self.pvaldf = None
        self.adjpvaldf = None
        self.cellgenedf = None
        self.sigCelltypedf = None

    def heatmapdf_create(self):
        '''
        This function will generate df for HeatMap by scaling the occurrencedf
        '''
        occurrencedf = self.occurrencedf
        transposedf = pd.DataFrame.transpose(occurrencedf)
        scaled_df = scale_dataframe(transposedf)
        scaled_df.columns = occurrencedf.index
        scaled_df = scaled_df.set_index(occurrencedf.columns)
        self.heatmapdf = scaled_df
        self.filter_heatmapdf()

    def filter_heatmapdf(self):
        '''
        This method filters rows and columns with sum < 1 in HeatMapdf
        '''
        df = self.heatmapdf
        self.filheatmapdf = df.loc[df.sum(1) > 1, df.sum(0) > 1]

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
        hclustHeatmap.path = os.path.sep.join([path, 'GCAM_heatMap.png'])
        fig, axm, axcb, cb = hclustHeatmap.plot()


    def fisher_occurrence_test(self):
        '''
        This method will calculate significance of celltypes per gene using their occurrence.
        Statistical test used is Fisher Exact Test
        '''
        occu_df = self.occurrencedf
        pvaldf = pd.DataFrame()#occu_df
        adjpvaldf = pd.DataFrame()#occu_df
        matsum = occu_df.sum().sum()
        for k, v in occu_df.iterrows():
            key = v.keys()
            #print v.shape
            rowsum = v.sum()
            for i in range(0, v.shape[0]):
                value = v[i]
                colsum = occu_df[[i]].sum()[0] - value
                rsum = rowsum - value
                if value != 0:
                    oddsratio, pval = stats.fisher_exact([[value, rsum], [colsum, matsum-(value+rsum+colsum)]])
                else:
                    pval = 1
                pvaldf.loc[k, key[i]] = pval
                if pval < 0.05/v.shape[0]:
                    adjpvaldf.loc[k, key[i]] = pval*v.shape[0]
                else:
                    adjpvaldf.loc[k, key[i]] = 1
        self.pvaldf = pvaldf
        self.adjpvaldf = adjpvaldf
        self.celltype_overrepresntation_list()  #### def()

    def celltype_overrepresntation_list(self):
        '''
        This method will save the result of significance in one DF.
        '''
        significance = 1
        column = ['CellType', 'Genes', 'P-val', 'FDR']
        cellgenedf = pd.DataFrame()
        for celltype, v in self.pvaldf.iterrows():
            for gene, pval in v.iteritems():
                if pval < significance:
                    cellgenedf = cellgenedf.append(pd.Series([celltype, gene, pval, self.adjpvaldf.loc[celltype, gene]])
                                                   , ignore_index=True)
        cellgenedf.columns = column
        self.cellgenedf = cellgenedf
        self.fisher_significant_celltypes()

    def fisher_significant_celltypes(self):
        '''
        This method will test the combined significance of celltype in the data and help predicts
        its association with user given data.
        '''
        sigcelltype = pd.DataFrame()
        #print self.cellgenedf
        cellgroup = self.cellgenedf.groupby(self.cellgenedf['CellType'])
        cellgenedf = self.cellgenedf
        c = len(cellgenedf[cellgenedf['FDR'] < 0.05])
        d = len(cellgenedf) #[cellgenedf['P-val'] < 0.5]
        for celltype, val in cellgroup:
            #print celltype
            #print val
            a = len(val[val['FDR'] < 0.05])
            b = len(val[val['P-val'] < 0.05]) - a
            cc = c - a
            dd = d - (a+b+cc)
            #print a, ':', b, ':', cc, ':', dd, c, d
            if a > 0:
                oddsRatio, p = stats.fisher_exact([[a, b], [cc, dd]])
                #print 'celltype:'+celltype, p
                sigcelltype = sigcelltype.append(pd.Series([celltype, a, p]), ignore_index=True)
        sigcelltype.columns = ['celltype', 'genecluster', 'P-val']
        length = len(sigcelltype)
        for k, v in sigcelltype.iterrows():
            if v['P-val'] < 0.05/length:
                sigcelltype.loc[k, 'FDR'] = v['P-val']*length
            else:sigcelltype.loc[k, 'FDR'] = 1
        self.sigCelltypedf = sigcelltype


def scale(val, src, dst):
    '''
    This returns scaled value
    :param val: value to be scaled
    :param src: max and min of values to be scaled
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
        rows = []
        for val in v:
            #print val
            old_min = min(v)
            old_max = max(v)
            if not isinstance(val, str):
                #print val
                rows.append(scale(val, (old_min, old_max), (new_min, new_max)))
        list_of_rows.append(rows)
    return pd.DataFrame(data=list_of_rows)