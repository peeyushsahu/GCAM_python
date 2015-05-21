from PubmedSearch.HclustHeatMap import HiearchicalHeatmap
from matplotlib import pyplot
__author__ = 'peeyush'
import pandas as pd
import scipy.stats as stats


class SignificanceObject():
    '''
    This Object will hold three dataframes
    '''
    def __init__(self, occurrencedf, heatmapdf=None, significancedf=None):
        self.occurrencedf = occurrencedf
        self.heatmapdf = heatmapdf
        self.significancedf = significancedf
        self.filheatmapdf = None
        self.pvaldf = None
        self.adjpvaldf = None

    def heatmap_create(self):
        occurrencedf = self.occurrencedf
        transposedf = pd.DataFrame.transpose(occurrencedf)
        scaled_df = scale_dataframe(transposedf)
        scaled_df.columns = occurrencedf.index
        scaled_df = scaled_df.set_index(occurrencedf.columns)
        self.heatmapdf = scaled_df
        self.filter_heatmapdf()

    def filter_heatmapdf(self):
        '''
        Filter rows and columns with sum < 1
        '''
        df = self.heatmapdf
        self.filheatmapdf = df.loc[df.sum(1) > 1, df.sum(0) > 1]

    def plot_heatmap(self):
        graph = TestHeatmap()
        graph.plotmap(self.filheatmapdf)
        return graph

    def fisher_test(self):
        occu_df = self.occurrencedf
        pvaldf = pd.DataFrame(occu_df)
        adjpvaldf = pd.DataFrame(occu_df)
        matsum = occu_df.sum().sum()
        for k, v in occu_df.iterrows():
            key = v.keys()
            #print v.shape
            rowsum = v.sum()
            for i in range(0, v.shape[0]):
                value = v[i]
                colsum = occu_df[[i]].sum()[0] - value
                rsum = rowsum - value
                oddsratio, pval = stats.fisher_exact([[value, rsum], [colsum, matsum-(value+rsum+colsum)]])
                pvaldf.loc[k, key[i]] = pval
                if pval < 0.05/v.shape[0]:
                    adjpvaldf.loc[k, key[i]] = pval*v.shape[0]
                else:
                    adjpvaldf.loc[k, key[i]] = 1
        self.pvaldf = pvaldf
        self.adjpvaldf = adjpvaldf

class TestHeatmap(HiearchicalHeatmap): ## This should create a object of HiearchicalHeatmap
    short_name = 'GCAM_HeatMap'

    def plotmap(self, df):
        self.frame = df
        self.path = '/home/peeyush/Desktop/' + self. short_name + '.pdf'
        fig, axm, axcb, cb = HiearchicalHeatmap.plot(self)
        cb.set_label("Random value")
        pyplot.savefig(self.path)


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
    print 'Process: scaling of dataframe'
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