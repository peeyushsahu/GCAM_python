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
        #print transposedf.head()
        scaled_df = scale_dataframe(transposedf)
        #print scaled_df.head()
        #print occurrencedf.index
        #print scaled_df
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
        hclustHeatmap.path = os.path.sep.join([path, 'GCAM_heatMap.svg'])
        fig, axm, axcb, cb = hclustHeatmap.plot()


    def fisher_occurrence_test(self):
        '''
        This method will calculate significance of celltypes per gene using their occurrence.
        Statistical test used is Fisher Exact Test
        '''
        occu_df = self.occurrencedf
        pvaldf = pd.DataFrame()
        adjpvaldf = pd.DataFrame()
        enrichmentdf = pd.DataFrame()
        matsum = occu_df.sum().sum()
        #print occu_df.head()
        for k, v in occu_df.iterrows():
            key = v.keys()
            rowsum = v.sum()
            for i in range(0, v.shape[0]):
                value = v[i]
                #print 'colsum',occu_df[[i]].sum()[0]
                if not value == 0:
                    colsum = occu_df[[i]].sum()[0] - value
                    rsum = rowsum - value
                    #print rowsum, value, colsum
                    #enrichment = scale(float(value), [min(v), max(v)], [0, 100])
                    enrichment = float(value)/occu_df[[i]].sum()[0]
                    if value != 0:
                        oddsratio, pval = stats.fisher_exact([[value, rsum], [colsum, matsum-(value+rsum+colsum)]])
                    else:
                        pval = 1
                    pvaldf.loc[k, key[i]] = pval
                    enrichmentdf.loc[k, key[i]] = enrichment
                    if pval < 0.05/v.shape[0]:
                        adjpvaldf.loc[k, key[i]] = pval*v.shape[0]
                    else:
                        adjpvaldf.loc[k, key[i]] = 1
                else:
                    pvaldf.loc[k, key[i]] = 1
                    adjpvaldf.loc[k, key[i]] = 1
                    enrichmentdf.loc[k, key[i]] = 0
        self.pvaldf = pvaldf
        self.adjpvaldf = adjpvaldf
        self.celltype_overrepresntation_list(enrichmentdf)  #### def()

    def celltype_overrepresntation_list(self, enrichmentdf):
        '''
        This method will save the result of significance in one DF.
        '''
        significance = 1
        column = ['celltype', 'gene', 'enrichment', 'p-val', 'q-val']
        cellgenedf = pd.DataFrame()
        for celltype, v in self.pvaldf.iterrows():
            for gene, pval in v.iteritems():
                if pval < significance:
                    cellgenedf = cellgenedf.append(pd.Series([celltype, gene, enrichmentdf.loc[celltype, gene], pval,
                                                         self.adjpvaldf.loc[celltype, gene]]), ignore_index=True)
        #print cellgenedf.head(10)
        cellgenedf.columns = column
        self.cellgenedf = cellgenedf
        #self.filter_cellgenedf()  # Filter single cell multigene enrihment
        self.fisher_significant_celltypes()

    def filter_cellgenedf(self):
        '''
        This method will remove one gene to multi celltype.
        '''
        group_cellgene = self.cellgenedf.groupby('gene')
        new_cellgenedf = pd.DataFrame(columns=self.cellgenedf.columns)
        siggenes = group_cellgene.groups
        for gene in siggenes:
            df = group_cellgene.get_group(gene).sort('enrichment', ascending=False)
            for ind, row in df.iterrows():
                #print ind, row
                if row['enrichment'] > 0.05 or row['p-val'] < 0.05:
                    new_cellgenedf = new_cellgenedf.append(row)
            '''
            if len(df) > 2:
                new_cellgenedf = new_cellgenedf.append(df.iloc[0, ])
                new_cellgenedf = new_cellgenedf.append(df.iloc[1, ])
                new_cellgenedf = new_cellgenedf.append(df.iloc[2, ])
            else:
                new_cellgenedf = new_cellgenedf.append(df)
            '''
        #print new_cellgenedf
        self.cellgenedf = new_cellgenedf


    def fisher_significant_celltypes(self):
        '''
        This method will test the combined significance of celltype in the data and help predicts
        its association with user given data.
        '''
        sigcelltype = pd.DataFrame()
        #print self.cellgenedf
        cellgroup = self.cellgenedf.groupby(self.cellgenedf['celltype'])
        cellgenedf = self.cellgenedf
        c = len(cellgenedf[cellgenedf['enrichment'] > 0.05])
        d = len(cellgenedf) #[cellgenedf['P-val'] < 0.5]
        for celltype, val in cellgroup:
            #print celltype
            if len(val[val['enrichment'] >= 0.1]) > 1:
                #print val
                a = len(val[val['enrichment'] >= 0.1])
                b = len(val) - a
                cc = c - a
                dd = d - (a+b+cc)
                #print a, ':', b, ':', cc, ':', dd, c, d
                if a > 0:
                    oddsRatio, p = stats.fisher_exact([[a, b], [cc, dd]])
                    #print 'celltype:'+celltype, a, p
                    sigcelltype = sigcelltype.append(pd.Series([celltype, a, p]), ignore_index=True)
        sigcelltype.columns = ['celltype', 'genecluster', 'p-val']
        length = len(sigcelltype[sigcelltype['p-val'] < 0.05])
        # Significant celltype check
        if length < 1:
            raise ValueError('No siginificant celltypes.')

        for k, v in sigcelltype.iterrows():
            if v['p-val'] < 0.05/length:
                sigcelltype.loc[k, 'q-val'] = v['p-val']*length
            else:sigcelltype.loc[k, 'q-val'] = 1
        self.sigCelltypedf = sigcelltype


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
