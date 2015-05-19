__author__ = 'peeyush'
import pandas as pd


class SignificanceObject():
    '''
    This Object will hold three dataframes
    '''
    def __init__(self, occurrencedf, heatmapdf=None, significancedf=None):
        self.occurrencedf = occurrencedf
        self.heatmapdf = heatmapdf
        self.significancedf = significancedf

    def heatmap_create(self):
        occurrencedf = self.occurrencedf
        transposedf = pd.DataFrame.transpose(occurrencedf)
        scaled_df = scale_dataframe(transposedf)
        scaled_df.columns = occurrencedf.index
        scaled_df = scaled_df.set_index(occurrencedf.columns)
        self.heatmapdf = scaled_df

    def plot_heatmap(df):
        df = df #self.heatmapdf
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap
        # Plot it out
        fig, ax = plt.subplots()
        custom_map = LinearSegmentedColormap.from_list(name='custom_div_cmap', colors =['blue', 'white', 'red'], N=299)
        heatmap = ax.pcolor(df, cmap=custom_map, alpha=0.5)
        # Format
        fig = plt.gcf()
        fig.set_size_inches(15, 11)
        # turn off the frame
        ax.set_frame_on(False)

        # put the major ticks at the middle of each cell
        #ax.set_yticks(np.arange(df.shape[0]) + 0.5, minor=False)
        #ax.set_xticks(np.arange(df.shape[1]) + 0.5, minor=False)

        # want a more natural, table-like display
        ax.invert_yaxis()
        ax.xaxis.tick_top()

        # note I could have used nba_sort.columns but made "labels" instead
        ax.set_xticklabels(df.columns, minor=False)
        ax.set_yticklabels(df.index, minor=False)
        # rotate the
        plt.xticks(rotation=90)

        ax.grid(False)

        # Turn off all the ticks
        ax = plt.gca()
        plt.savefig('/home/peeyush/Desktop/new_heatmap_p6.png')

def hclustering(df):
    '''
    
    :param df:
    :return:
    '''
    import fastcluster
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import dendrogram
    from scipy.spatial.distance import pdist
    import scipy.cluster.hierarchy as hierarchy
    Z = fastcluster.linkage(X, method ='single')
    leaves = hierarchy.leaves_list(Z)
    count = 0
    for k, v in df.iterrows():
        df.loc[k, 'cluster'] = leaves[count]
        count += 1


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