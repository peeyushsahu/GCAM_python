__author__ = 'peeyush'


"""
This code was adapted from the following recipe:
    * http://altanalyze.blogspot.se/2012/06/hierarchical-clustering-heatmaps-in.html
    * http://code.activestate.com/recipes/578175/

Which was in turn inspired by many other posts:
   * http://stackoverflow.com/questions/7664826
   * http://stackoverflow.com/questions/2982929
   * http://stackoverflow.com/questions/2455761

Running this with cosine or other distance metrics can often produce negative Z scores during clustering, so adjustments to the clustering may be required. Information about distance measures can be found here:
   * http://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html
   * http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.cdist.html

The documentation about the custom color gradients can be found here:
   * http://matplotlib.sourceforge.net/examples/pylab_examples/custom_cmap.html
"""

# Third party modules #
import numpy, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dist
import os
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd


###############################################################################
# Create Custom Color Gradients #
red_black_sky = {'red':   ((0.0, 0.0, 0.0), (0.5, 0.0, 0.1), (1.0, 1.0, 1.0)),
                     'green': ((0.0, 0.0, 0.9), (0.5, 0.1, 0.0), (1.0, 0.0, 0.0)),
                     'blue':  ((0.0, 0.0, 1.0), (0.5, 0.1, 0.0), (1.0, 0.0, 0.0))}
red_black_blue = {'red':   ((0.0, 0.0, 0.0), (0.5, 0.0, 0.1), (1.0, 1.0, 1.0)),
                     'green': ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)),
                     'blue':  ((0.0, 0.0, 1.0), (0.5, 0.1, 0.0), (1.0, 0.0, 0.0))}
red_black_green = {'red':   ((0.0, 0.0, 0.0), (0.5, 0.0, 0.1), (1.0, 1.0, 1.0)),
                     'blue':  ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)),
                     'green': ((0.0, 0.0, 1.0), (0.5, 0.1, 0.0), (1.0, 0.0, 0.0))}
yellow_black_blue = {'red':   ((0.0, 0.0, 0.0), (0.5, 0.0, 0.1), (1.0, 1.0, 1.0)),
                     'green': ((0.0, 0.0, 0.8), (0.5, 0.1, 0.0), (1.0, 1.0, 1.0)),
                     'blue':  ((0.0, 0.0, 1.0), (0.5, 0.1, 0.0), (1.0, 0.0, 0.0))}

make_cmap = lambda x: matplotlib.colors.LinearSegmentedColormap('my_colormap', x, 256)
color_gradients = {'red_black_sky'      : make_cmap(red_black_sky),
                   'red_black_blue'     : make_cmap(red_black_blue),
                   'red_black_green'    : make_cmap(red_black_green),
                   'yellow_black_blue'  : make_cmap(yellow_black_blue),
                   'red_white_blue'     : pyplot.cm.bwr,
                   'seismic'            : pyplot.cm.seismic,
                   'green_white_purple' : pyplot.cm.PiYG_r,
                   'coolwarm'           : pyplot.cm.coolwarm,}

###############################################################################
class HiearchicalHeatmap():
    """A common use case for biologists analyzing their gene expression data is to cluster and visualize patterns of expression in the form of a heatmap and associated dendrogram."""
    def __init__(self):
        self.row_method = 'single'     # Can be: linkage, single, complete, average, weighted, centroid, median, ward
        self.column_method = 'single'     # Can be: linkage, single, complete, average, weighted, centroid, median, ward
        self.row_metric = 'braycurtis' # Can be: see scipy documentation
        self.column_metric = 'braycurtis' # Can be: see scipy documentation
        self.gradient_span = 'only_max'   # Can be: min_to_max, min_to_max_centered, only_max, only_min
        self.color_gradient = 'red_white_blue'   # Can be: see color_gradients dictionary
        self.fig_weight = 12
        self.fig_height = 8.5
        self.frame = None
        self.path = None

    def plot(self):
        # Names #
        row_header = self.frame.index
        column_header = self.frame.columns

        # What color to use #
        cmap = color_gradients[self.color_gradient]

        # Scale the max and min colors #
        value_min = self.frame.min().min()
        value_max = self.frame.max().max()
        if self.gradient_span == 'min_to_max_centered':
            value_max = max([value_max, abs(value_min)])
            value_min = value_max * -1
        if self.gradient_span == 'only_max': value_min = 0
        if self.gradient_span == 'only_min': value_max = 0
        norm = matplotlib.colors.Normalize(value_min, value_max)

        # Scale the figure window size #
        fig = pyplot.figure(figsize=(self.fig_weight, self.fig_height))

        # Calculate positions for all elements #
        # ax1, placement of dendrogram 1, on the left of the heatmap
        ### The second value controls the position of the matrix relative to the bottom of the view
        [ax1_x, ax1_y, ax1_w, ax1_h] = [0.05, 0.22, 0.2, 0.6]
        width_between_ax1_axr = 0.004
        ### distance between the top color bar axis and the matrix
        height_between_ax1_axc = 0.004
        ### Sufficient size to show
        color_bar_w = 0.015

        # axr, placement of row side colorbar #
        ### second to last controls the width of the side color bar - 0.015 when showing
        [axr_x, axr_y, axr_w, axr_h] = [0.31, 0.1, color_bar_w, 0.6]
        axr_x = ax1_x + ax1_w + width_between_ax1_axr
        axr_y = ax1_y; axr_h = ax1_h
        width_between_axr_axm = 0.004

        # axc, placement of column side colorbar #
        ### last one controls the hight of the top color bar - 0.015 when showing
        [axc_x, axc_y, axc_w, axc_h] = [0.4, 0.63, 0.5, color_bar_w]
        axc_x = axr_x + axr_w + width_between_axr_axm
        axc_y = ax1_y + ax1_h + height_between_ax1_axc
        height_between_axc_ax2 = 0.004

        # axm, placement of heatmap for the data matrix #
        [axm_x, axm_y, axm_w, axm_h] = [0.4, 0.9, 2.5, 0.5]
        axm_x = axr_x + axr_w + width_between_axr_axm
        axm_y = ax1_y; axm_h = ax1_h
        axm_w = axc_w

        # ax2, placement of dendrogram 2, on the top of the heatmap #
        ### last one controls hight of the dendrogram
        [ax2_x, ax2_y, ax2_w, ax2_h] = [0.3, 0.72, 0.6, 0.15]
        ax2_x = axr_x + axr_w + width_between_axr_axm
        ax2_y = ax1_y + ax1_h + height_between_ax1_axc + axc_h + height_between_axc_ax2
        ax2_w = axc_w

        # axcb - placement of the color legend #
        [axcb_x, axcb_y, axcb_w, axcb_h] = [0.07, 0.88, 0.18, 0.09]

        # Compute and plot top dendrogram #
        if self.column_method:
            d2 = dist.pdist(self.frame.transpose())
            D2 = dist.squareform(d2)
            ax2 = fig.add_axes([ax2_x, ax2_y, ax2_w, ax2_h], frame_on=True)
            Y2 = sch.linkage(D2, method=self.column_method, metric=self.column_metric)
            Z2 = sch.dendrogram(Y2)
            ind2 = sch.fcluster(Y2, 0.4*max(Y2[:,2]), 'distance')
            ax2.set_xticks([])
            ax2.set_yticks([])
            ### apply the clustering for the array-dendrograms to the actual matrix data
            idx2 = Z2['leaves']
            self.frame = self.frame.iloc[:,idx2]
            ### reorder the flat cluster to match the order of the leaves the dendrogram
            ind2 = ind2[idx2]
        else: idx2 = range(self.frame.shape[1])

        # Compute and plot left dendrogram #
        if self.row_method:
            d1 = dist.pdist(self.frame)
            D1 = dist.squareform(d1)
            ax1 = fig.add_axes([ax1_x, ax1_y, ax1_w, ax1_h], frame_on=True)
            Y1 = sch.linkage(D1, method=self.row_method, metric=self.row_metric)
            Z1 = sch.dendrogram(Y1, orientation='right')
            ind1 = sch.fcluster(Y1, 0.4*max(Y1[:,2]), 'distance')
            ax1.set_xticks([])
            ax1.set_yticks([])
            ### apply the clustering for the array-dendrograms to the actual matrix data
            idx1 = Z1['leaves']
            self.frame = self.frame.iloc[idx1,:]
            ### reorder the flat cluster to match the order of the leaves the dendrogram
            ind1 = ind1[idx1]
        else: idx1 = range(self.frame.shape[0])

        # Plot distance matrix #
        axm = fig.add_axes([axm_x, axm_y, axm_w, axm_h])
        axm.matshow(self.frame, aspect='auto', origin='lower', cmap=cmap, norm=norm)
        axm.set_xticks([])
        axm.set_yticks([])

        # Add text #
        new_row_header = []
        new_column_header = []
        for i in range(self.frame.shape[0]):
            axm.text(self.frame.shape[1]-0.5, i, '  ' + row_header[idx1[i]], verticalalignment="center")
            new_row_header.append(row_header[idx1[i]] if self.row_method else row_header[i])
        for i in range(self.frame.shape[1]):
            axm.text(i, -0.55, ' '+column_header[idx2[i]], rotation=90, verticalalignment="top", horizontalalignment="center")
            new_column_header.append(column_header[idx2[i]] if self.column_method else column_header[i])

        # Plot column side colorbar #
        if self.column_method:
            axc = fig.add_axes([axc_x, axc_y, axc_w, axc_h])
            cmap_c = matplotlib.colors.ListedColormap(['r', 'g', 'b', 'y', 'w', 'k', 'm'])
            dc = numpy.array(ind2, dtype=int)
            dc.shape = (1,len(ind2))
            axc.matshow(dc, aspect='auto', origin='lower', cmap=cmap_c)
            axc.set_xticks([])
            axc.set_yticks([])

        # Plot column side colorbar #
        if self.row_method:
            axr = fig.add_axes([axr_x, axr_y, axr_w, axr_h])
            dr = numpy.array(ind1, dtype=int)
            dr.shape = (len(ind1),1)
            cmap_r = matplotlib.colors.ListedColormap(['r', 'g', 'b', 'y', 'w', 'k', 'm'])
            axr.matshow(dr, aspect='auto', origin='lower', cmap=cmap_r)
            axr.set_xticks([])
            axr.set_yticks([])

        # Plot color legend #
        ### axes for colorbar
        axcb = fig.add_axes([axcb_x, axcb_y, axcb_w, axcb_h], frame_on=False)
        cb = matplotlib.colorbar.ColorbarBase(axcb, cmap=cmap, norm=norm, orientation='horizontal')
        axcb.set_title("colorkey", fontsize=4)
        max_cb_ticks = 5
        axcb.xaxis.set_major_locator(pyplot.MaxNLocator(max_cb_ticks))

        # Render the graphic #
        if len(row_header) > 30 or len(column_header) > 30:
            #print ('more than 80', len(row_header))
            pyplot.rcParams['font.size'] = 2
        else:
            #print 'less than 80', len(row_header)
            pyplot.rcParams['font.size'] = 6
        #print(pyplot.rcParams.find_all('\.size'))
        cb.set_label("Enrichment scale", fontsize=8)
        pyplot.savefig(self.path)
        pyplot.clf()
        # Return figure #
        return fig, axm, axcb, cb


def radar_factory(num_vars, frame='circle'):
    import numpy as np
    #import seaborn as sns
    #sns.set_style("white")
    import matplotlib.pyplot as plt
    from matplotlib.path import Path
    from matplotlib.spines import Spine
    from matplotlib.projections.polar import PolarAxes
    from matplotlib.projections import register_projection
    """Create a radar chart with `num_vars` axes.

    This function creates a RadarAxes projection and registers it.

    Parameters
    ----------
    num_vars : int
        Number of variables for radar chart.
    frame : {'circle' | 'polygon'}
        Shape of frame surrounding axes.

    """
    # calculate evenly-spaced axis angles
    theta = np.linspace(0, 2*np.pi, num_vars, endpoint=False)
    # rotate theta such that the first axis is at the top
    theta += np.pi/2

    def draw_poly_patch(self):
        verts = unit_poly_verts(theta)
        return plt.Polygon(verts, closed=True, edgecolor='k')

    def draw_circle_patch(self):
        # unit circle centered on (0.5, 0.5)
        return plt.Circle((0.5, 0.5), 0.5)

    patch_dict = {'polygon': draw_poly_patch, 'circle': draw_circle_patch}
    if frame not in patch_dict:
        raise ValueError('unknown value for `frame`: %s' % frame)

    class RadarAxes(PolarAxes):

        name = 'radar'
        # use 1 line segment to connect specified points
        RESOLUTION = 1
        # define draw_frame method
        draw_patch = patch_dict[frame]

        def fill(self, *args, **kwargs):
            """Override fill so that line is closed by default"""
            closed = kwargs.pop('closed', True)
            return super(RadarAxes, self).fill(closed=closed, *args, **kwargs)

        def plot(self, *args, **kwargs):
            """Override plot so that line is closed by default"""
            lines = super(RadarAxes, self).plot(*args, **kwargs)
            for line in lines:
                self._close_line(line)

        def _close_line(self, line):
            x, y = line.get_data()
            # FIXME: markers at x[0], y[0] get doubled-up
            if x[0] != x[-1]:
                x = np.concatenate((x, [x[0]]))
                y = np.concatenate((y, [y[0]]))
                line.set_data(x, y)

        def set_varlabels(self, labels):
            self.set_thetagrids(np.degrees(theta), labels)

        def _gen_axes_patch(self):
            return self.draw_patch()

        def _gen_axes_spines(self):
            if frame == 'circle':
                return PolarAxes._gen_axes_spines(self)
            # The following is a hack to get the spines (i.e. the axes frame)
            # to draw correctly for a polygon frame.

            # spine_type must be 'left', 'right', 'top', 'bottom', or `circle`.
            spine_type = 'circle'
            verts = unit_poly_verts(theta)
            # close off polygon by repeating first vertex
            verts.append(verts[0])
            path = Path(verts)

            spine = Spine(self, spine_type, path)
            spine.set_transform(self.transAxes)
            return {'polar': spine}

    register_projection(RadarAxes)
    return theta


def unit_poly_verts(theta):
    """Return vertices of polygon for subplot axes.

    This polygon is circumscribed by a unit circle centered at (0.5, 0.5)
    """
    x0, y0, r = [0.5] * 3
    verts = [(r*np.cos(t) + x0, r*np.sin(t) + y0) for t in theta]
    return verts


def plot_radar(gene2plot, tissue, path):
    '''
    Plot radar plot
    :param data:
    :param path:
    :return:
    '''
    data = gene2plot
    #[0.00364299,0.0582878,0.04189435,0.13661202,0.10928962,0.14754098,0.00728597,0.0582878,0.40801457,0.3]
    spoke_labels = tissue
    #['lung', 'bone', 'blood', 'epithelial', 'progenitor', 'nervous', 'endocrine', 'epidermal', 'immune','liver']
    N = len(data)
    theta = radar_factory(N, frame='polygon')

    plt.rcParams['font.size'] = 10
    fig = plt.figure(figsize=(6, 6))
    fig.subplots_adjust(wspace=0.25, hspace=0.20, top=0.85, bottom=0.05)

    ax = fig.add_subplot(1, 1, 1, projection='radar')
    #plt.rgrids([0.1, 0.2, 0.3, 0.4])
    #ax.set_title('Fraction of gene', weight='bold', fontsize=15, position=(0.5, 1.1),
    #             horizontalalignment='center', verticalalignment='center')
    ax.plot(theta, data, color='r', lw=1)
    ax.fill(theta, data, facecolor='r', alpha=0.25)
    ax.set_varlabels(spoke_labels)
    plt.figtext(0.5, 0.965, 'Gene enrichment fraction for (tissue)* types',
                ha='center', color='black', weight='bold', fontsize=15)
    plt.savefig(os.path.join(path, 'GCAM_redar.svg'))
    plt.clf()
    plt.close()
    #plt.show()


def stack_barplot(sigCelltypedf, args, plotdf, path):
    import seaborn as sns
    sns.set_context("paper")
    sns.despine()
    plotDF = pd.DataFrame()
    if not args.key_celltype_list:
        sigCelltypedf = sigCelltypedf.sort_values('p-val', ascending=True)
        celltypes = 0
        for k, v in sigCelltypedf.iterrows():
            if (v['genecluster'] > 20) and (celltypes < 15) and (v['celltype'] in list(plotdf.index)):
                plotDF = plotDF.append(plotdf.loc[v['celltype']])
            celltypes += 1
    else:
        plotDF = plotdf
    #print(plotDF)
    plotDF = plotDF.sort_index(axis=0, ascending=True)
    # Join three palette in case of more than 8 items
    if len(plotDF) > 8:
        d_colors = sns.hls_palette(8, l=.3, s=.7)
        d_colors.pop(3)
        d_colors1 = sns.hls_palette(8, l=.5, s=1)
        d_colors1.pop(3)
        d_colors2 = sns.hls_palette(8, l=.7, s=.8)
        d_colors2.pop(3)
        n_d_colors = []
        for i in range(0,len(d_colors)):
            n_d_colors.append(d_colors[i])
            n_d_colors.append(d_colors1[i])
            n_d_colors.append(d_colors2[i])
            #print(len(n_d_colors))
    else:
        n_d_colors = sns.color_palette("hls", 8)

    # Scaling plotDF for stack bar plotting to 1
    #print(plotDF)
    df1 = plotDF.div(plotDF.sum(axis=0), axis=1)
    # Initialize the vertical-offset for the stacked bar chart.
    y_offset = np.array([0.0] * len(df1.columns))
    ind = np.arange(len(df1.columns))    # the x locations for the groups
    width = 0.5 # the width of the bars: can also be len(x) sequence
    plt.figure(figsize=[10,6])
    handls = []
    i = 0
    for k, v in df1.iterrows():
        h = plt.bar(ind+0.05, list(v), width, bottom=y_offset, color=n_d_colors[i], edgecolor='w')#alpha=0.7
        handls.append(h[0])
        y_offset = y_offset + v
        i += 1
    plt.ylabel('Cell fractions')
    plt.xlabel('Samples')
    plt.title("Relative change in celltype expression")
    plt.ylim(0,1)
    plt.xlim(-0.2,len(df1.columns))
    plt.xticks(ind + width/2. + 0.05, list(df1.columns), rotation=45)
    plt.yticks(np.arange(0, 1, 0.1))
    lgd = plt.legend(tuple(handls[::-1]), tuple(list(df1.index)), loc='center right', bbox_to_anchor=(1.2, 0.5))
    plt.tick_params(direction='out')
    plt.tight_layout()
    plt.savefig(os.path.join(path, 'GCAM_stacks.svg'), bbox_extra_artists=(lgd, ), bbox_inches='tight')
    plt.clf()
    plt.close()


def heatmap_Sigcelltype(args, df, path):
    import seaborn as sns
    '''
    To plot stack plot df as heatmap
    '''
    #print(args.key_celltype_list)
    if args.key_celltype_list:
        cell_type = ['macrophage', 'Alveolar macrophage',
             'm1 macrophage','m2 macrophage', 'monocyte', 'dendritic cell', 'glial cell',
             'neutrophil', 'mast cell', 'Natural killer cell', 'Kupffer cell', 'Plasma cell',
             'eosinophil', 'naive B cell', 'memory B cell', 'B lymphocyte', 'T lymphocyte',
             'naive T cell', 'memory T cell', 'CD8 T cell', 'CD4 T cell', 'regulatory T cell','Cytotoxic T cell',
             'helper T cell']
        # creating df for heatmap
        new_df = pd.DataFrame(0, columns=df.columns, index=cell_type)
        #print(new_df)
        for k, v in df.iterrows():
            for c, val in v.iteritems():
                #print(c, val)
                new_df.loc[k, c] = val
        # plotting df
        new_df = new_df.T
        sns.set_context("talk")
        cmap = sns.diverging_palette(255, 15, sep=20, n=3, as_cmap=True)
        plt.clf()
        plt.figure(figsize=[20,10])
        sns.heatmap(new_df.round(2), cmap = cmap, vmin=0, vmax=0.2, yticklabels=True, cbar=False,
            xticklabels=True, linecolor='#ffffff',linewidths=0.01, square=True, annot=True)
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(path, 'GCAM_cofficients.pdf'))
        plt.close()
    else:
        # creating df for heatmap
        df = df.T
        sns.set_context("talk")
        cmap = sns.diverging_palette(255, 15, sep=20, n=3, as_cmap=True)
        plt.clf()
        plt.figure(figsize=[20,10])
        sns.heatmap(df.round(2), cmap = cmap, vmin=0, vmax=0.2, yticklabels=True, cbar=False,
            xticklabels=True, linecolor='#ffffff',linewidths=0.01, square=True, annot=True)
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(path, 'GCAM_cofficients.pdf'))
        plt.close()

#########################################


def plot_celltypesignificance(path, plotdf, args):
    import matplotlib.pyplot as plt
    import numpy as np
    import sys
    print ('plotting celltype significance plot')
    plotdf = plotdf[(plotdf['genecluster'] >= 5)]
    if len(plotdf) < 1:
        sys.exit('Not enough genes for significant celltype plot, please see results at the GCAM_sigenes.xls')
    if args.subcommand_name == 'exprbased':
        plotdf = plotdf.sort_values('p-val',ascending=True)
    l = np.log2(plotdf['genecluster'].tolist())
    t = plotdf['p-val'].tolist()
    s = range(1, len(plotdf)+1)
    name = plotdf['celltype'].tolist()
    area = [abs((np.log10(x)) * 15) for x in t]
    color = np.random.random(len(t))
    #print area
    #plt.figure()
    plt.scatter(s, l, s=area, c=color, alpha=0.5)
    plt.grid(True, linestyle=':', color='black')

    # draw a thick red hline at y=0 that spans the xrange
    thres = np.log2(args.celltypeClusterSize)
    if thres < min(l): thres = min(l)
    h = plt.axhline(linewidth=1, color='r', y=thres, linestyle='--') #y=args.celltypeClusterSize

    for i in range(0, len(t)):
        plt.annotate(name[i], xy=(s[i], l[i]), xycoords='data',
            xytext=(0,15), textcoords='offset points',
            ha='center', va='bottom',
            bbox=dict(boxstyle='round, pad=0.2', fc='yellow', alpha=0.2),
            fontsize=10)
## plot legend
    l1 = plt.scatter([],[], s=50, c='gray', alpha=0.5)
    l2 = plt.scatter([],[], s=200, c='gray', alpha=0.5)
    labels = ["less significant", "highly significant"]
    plt.legend([l1, l2], labels, ncol=2, frameon=True, fontsize=8,
    handlelength=2, loc = 4, borderpad = 0.5,
    handletextpad=1, scatterpoints = 1)

    plt.tick_params(axis='both', labelsize=8)
    plt.xlim(0, len(plotdf)+2)
    plt.ylim(min(l)-0.2, max(l)+0.2)
    plt.title('Cell-type significance', fontsize=14)
    plt.xlabel('Celltypes', fontsize=12)
    plt.ylabel('log2 Gene cluster size', fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(path, 'GCAM_SigCelltype.svg'))
    plt.clf()

'''
def heatmap_Sigcelltype(df, path, edgecolors='w', log=False):

    :param df:
    :param edgecolors:
    :param log:
    :return:
    import matplotlib.colors as mcolors
    import numpy as np
    import matplotlib.pyplot as plt

    width = len(df.columns)/7*10
    height = len(df.index)/7*10

    fig, ax = plt.subplots(figsize=(20,10))#(figsize=(width,height))
    dfMax = max(df.max()) + max(df.max())/15
    dfSize = np.linspace(0, dfMax, num=15)
    #print dfSize
    cmap, norm = mcolors.from_levels_and_colors(dfSize, ['#3e7d00', '#579619', '#a2c57f', '#b4d099', '#c7dcb2',
                                                         '#d9e7cc', '#ecf3e5', '#ffffff', '#fae5e5', '#f5cccc',
                                                         '#f0b2b2', '#eb9999', '#e67f7f', '#d21919'] ) # ['MidnightBlue', Teal]['Darkgreen', 'Darkred']
    heatmap = ax.pcolor(df,
                        edgecolors=edgecolors,  # put white lines between squares in heatmap
                        cmap=cmap,
                        norm=norm)
    data = df.values
    for y in range(data.shape[0]):
        for x in range(data.shape[1]):
            plt.text(x + 0.5 , y + 0.5, '%.4f' % data[y, x], #data[y,x] +0.05 , data[y,x] + 0.05
                 horizontalalignment='center',
                 verticalalignment='center',
                 size=10,
                 color='b')

    ax.autoscale(tight=True)  # get rid of whitespace in margins of heatmap
    ax.set_aspect('equal')  # ensure heatmap cells are square
    ax.xaxis.set_ticks_position('top')  # put column labels at the top
    ax.tick_params(bottom='off', top='off', left='off', right='off')  # turn off ticks

    ax.set_yticks(np.arange(len(df.index)) + 0.5)
    ax.set_yticklabels(df.index, size=15)
    ax.set_xticks(np.arange(len(df.columns)) + 0.5)
    ax.set_xticklabels(df.columns, rotation=90, size= 15)

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    #cax = divider.append_axes("right", "3%", pad="1%")
    #fig.colorbar(heatmap, cax=cax)
    plt.tight_layout()
    plt.savefig(os.path.join(path, 'GCAM_heatmap_SigCelltype.svg'))
    plt.clf()
    plt.close()
'''