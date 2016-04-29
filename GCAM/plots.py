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


def stack_barplot(args, sigcelltype, path, name='_', method='genebased'):
    #print 'Key_celltypes', key_celltypes
    d_labels = ["celltype"]
    d_widths = [0.1]
    stackplot_data = []
    d_colors = ['#006699', '#33D6AD', '#CCFF33', '#0033CC', '#FFFF00', '#009900', '#CC3300', '#FFCC99',
                '#99CCFF', '#CC6699', '#CC99FF', '#6600FF', '#CD853F', '#FFC0CB', '#FF9900', '#B0E0E6',
                '#800080', '#663399', '#70dd3e', '#BC8F8F', '#4169E1', '#8B4513', '#FA8072', '#F4A460',
                '#2E8B57', '#ffc8a2','#00FF00', '#FFFF66', '#191919', '#bcffad']
    if method == 'genebased' and not args.key_celltypes:
        if len(sigcelltype) > 15: sigcelltype = sigcelltype[:15]
        cells = [1.0]*len(sigcelltype)
        number_label = float(len(sigcelltype))
        adjustment = 0.02
        plot_color = d_colors[:len(cells)]
        #print sigcelltype
        #if method == 'genebased':
        genes = [0.0]*len(cells)
        stackplot_data.append(cells)
        sigcelltype = sigcelltype.sort(['P-val'], ascending=True)
        sigcelltype.index = range(0, len(sigcelltype))
        cell_type = sigcelltype['celltype']
        for ind, row in sigcelltype.iterrows():
            genes[ind] = float(row['genecluster'])
        stackplot_data.insert(0, genes)
        d_labels.insert(0, 'sample')
        d_widths.insert(0, 0.5)

    if not args.key_celltype_list and method == 'exprbased':
        #if len(sigcelltype) > 15: sigcelltype = sigcelltype[:15]
        cells = [1.0]*len(sigcelltype)
        number_label = float(len(sigcelltype))
        adjustment = 0.02
        plot_color = d_colors[:len(cells)]
        cell_type = list(sigcelltype.index)
        #print cell_type
        #print sigcelltype
        stackplot_data.append(cells)
        #cols = [col for col in sigcelltype.columns if 'FC' in col]
        #print sigcelltype.columns
        for col in sigcelltype.columns:
            genes = [0.0]*len(cells)
            for i, r in sigcelltype.iterrows():
                #print i
                ind = cell_type.index(i)
                genes[ind] = float(r[col])
            stackplot_data.insert(0, genes)
            d_labels.insert(0, col)
            d_widths.insert(0, 0.5)
            #print stackplot_data
    else:
        cell_type = ['macrophage', 'Alveolar macrophage',
                     'm1 macrophage','m2 macrophage', 'monocyte', 'dendritic cell', 'glial cell',
                     'neutrophil', 'mast cell', 'Natural killer cell', 'Kupffer cell', 'Plasma cell',
                     'eosinophil', 'naive B cell', 'memory B cell', 'B lymphocyte', 'T lymphocyte',
                     'naive T cell', 'memory T cell', 'CD8 T cell', 'CD4 T cell', 'regulatory T cell','Cytotoxic T cell',
                     'helper T cell']
        d_colors = ['#009688', '#35a79c', '#54b2a9', '#83d0c9', '#ffdaaf',
                    '#ff9237', '#eb5300', '#8d8282', '#bbbbbb', '#f5f8ee', '#feffc3',
                    '#fcff83', '#f9e909', '#fe2e2e', '#cb2424', '#b62020', '#e18f8f', '#fcbbbb',
                    '#efbbff','#d896ff', '#be29ec','#800080', '#660066', '#470047']

        cells = [1.0]*len(cell_type)
        number_label = float(len(cell_type))
        adjustment = 0.02
        plot_color = d_colors[:len(cells)]
        cell_type = list(map(str.lower, cell_type))

        if method == 'genebased':
            stackplot_data.append(cells)
            genes = [0.0]*len(cells)
            for i, r in sigcelltype.iterrows():
                ind = cell_type.index(r['celltype'])
                genes[ind] = float(r['genecluster'])
            stackplot_data.insert(0, genes)
            d_labels.insert(0, 'sample')
            d_widths.insert(0, 0.5)
        if method == 'exprbased':
            stackplot_data.append(cells)
            #cols = [col for col in sigcelltype.columns if 'FC' in col]
            for col in sigcelltype.columns:
                genes = [0.0]*len(cells)
                for i, r in sigcelltype.iterrows():
                    #print(i)
                    ind = cell_type.index(i)
                    genes[ind] = float(r[col])
                stackplot_data.insert(0, genes)
                d_labels.insert(0, col)
                d_widths.insert(0, 0.5)
    #print stackplot_data
    #print d_labels
    stplot = np.array(stackplot_data) #[genes, cells] stackplot list of list
    gap = 0.01
    fig = plt.figure(figsize=(11,8))
    ax6 = fig.add_subplot(111)
    stackedBarPlot(ax6,
                    stplot,
                    plot_color,
                    cell_type,
                    number_label,
                    adjustment,
                    edgeCols=['#ffffff']*int(number_label),
                    xLabels=d_labels,
                    scale=True,
                    gap=gap,
                    endGaps=True,
                    widths=d_widths
                    )
    plt.title("Relative abundance of celltype")
    plt.setp(ax6.xaxis.get_majorticklabels(), rotation=90)
    #fig.subplots_adjust(bottom=0.4)
    plt.tight_layout()
    plt.savefig(os.path.join(path, 'GCAM_'+method+'_stacks'+name+'.svg'))
    plt.clf()
    #return


def stackedBarPlot(ax,                                 # axes to plot onto
                   data,                               # data to plot
                   cols,                               # colors for each level
                   cell_type,
                   no_label,
                   adjustment,
                   xLabels = None,                     # bar specific labels
                   yTicks = 6.0,                        # information used for making y ticks ["none", <int> or [[tick_pos1, tick_pos2, ... ],[tick_label_1, tick_label2, ...]]
                   edgeCols=None,                      # colors for edges
                   showFirst=-1,                       # only plot the first <showFirst> bars
                   scale=False,                        # scale bars to same height
                   widths=None,                        # set widths for each bar
                   heights=None,                       # set heights for each bar
                   ylabel='',                          # label for x axis
                   xlabel='',                          # label for y axis
                   gap=0.,                             # gap between bars
                   endGaps=False                       # allow gaps at end of bar chart (only used if gaps != 0.)
                   ):

# data fixeratering

    # make sure this makes sense
    if showFirst != -1:
        showFirst = np.min([showFirst, np.shape(data)[0]])
        data_copy = np.copy(data[:showFirst]).transpose().astype('float')
        data_shape = np.shape(data_copy)
        if heights is not None:
            heights = heights[:showFirst]
        if widths is not None:
            widths = widths[:showFirst]
        showFirst = -1
    else:
        data_copy = np.copy(data).transpose()
    data_shape = np.shape(data_copy)

    # determine the number of bars and corresponding levels from the shape of the data
    num_bars = data_shape[1]
    levels = data_shape[0]

    if widths is None:
        widths = np.array([1] * num_bars)
        x = np.arange(num_bars)
    else:
        x = [0]
        for i in range(1, len(widths)):
            x.append(x[i-1] + (widths[i-1] + widths[i])/2)

    # stack the data --
    # replace the value in each level by the cumulative sum of all preceding levels
    data_stack = np.reshape([float(i) for i in np.ravel(np.cumsum(data_copy, axis=0))], data_shape)

    # scale the data is needed
    if scale:
        data_copy /= data_stack[levels-1]
        data_stack /= data_stack[levels-1]
        if heights is not None:
            print("WARNING: setting scale and heights does not make sense.")
            heights = None
    elif heights is not None:
        data_copy /= data_stack[levels-1]
        data_stack /= data_stack[levels-1]
        for i in np.arange(num_bars):
            data_copy[:,i] *= heights[i]
            data_stack[:,i] *= heights[i]

#------------------------------------------------------------------------------
# ticks

    if yTicks is not "none":
        # it is either a set of ticks or the number of auto ticks to make
        real_ticks = True
        try:
            k = len(yTicks[1])
        except:
            real_ticks = False

        if not real_ticks:
            yTicks = float(yTicks)
            if scale:
                # make the ticks line up to 100 %
                y_ticks_at = np.arange(yTicks)/(yTicks-1)
                y_tick_labels = np.array(["%0.2f"%(i * 100) for i in y_ticks_at])
            else:
                # space the ticks along the y axis
                y_ticks_at = np.arange(yTicks)/(yTicks-1)*np.max(data_stack)
                y_tick_labels = np.array([str(i) for i in y_ticks_at])
            yTicks=(y_ticks_at, y_tick_labels)

#------------------------------------------------------------------------------
# plot

    if edgeCols is None:
        edgeCols = ["none"]*len(cols)

    # take cae of gaps
    gapd_widths = [i - gap*2 if(widths.index(i) == len(widths)-1) else i-gap*2 for i in widths]
    #print gapd_widths

    # bars
    ax.bar(x,
           data_stack[0],
           color=cols[0],
           edgecolor=edgeCols[0],
           width=gapd_widths,
           linewidth=0.5,
           align='center'
           )

    for i in np.arange(1,levels):
        ax.bar(x,
               data_copy[i],
               bottom=data_stack[i-1],
               color=cols[i],
               edgecolor=edgeCols[i],
               width=gapd_widths,
               linewidth=0.5,
               align='center'
               )

    # borders
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    #ax.spines["left"].set_visible(False)

    # make ticks if necessary
    if yTicks is not "none":
        ax.tick_params(axis='y', which='both', labelsize=8, direction="out")
        ax.yaxis.tick_left()
        plt.yticks(yTicks[0], ['0%','20%','40%','60%','80%','100%'])

        # second y-axis for cell_type annotation
        ax2 = ax.twinx()
        ax2.tick_params(axis='y', which='both', labelsize=10, direction='in')
        ax2.yaxis.tick_right()
        plt.yticks(np.arange(no_label)/(no_label)+adjustment, cell_type)
    else:
        plt.yticks([], [])

    if xLabels is not None:
        ax.tick_params(axis='x', which='both', labelsize=8, direction="out")
        ax.xaxis.tick_bottom()
        plt.xticks(x, xLabels)
    else:
        plt.xticks([], [])

    # limits
    if endGaps:
        ax.set_xlim(-1.*widths[0]/2. - gap/2., np.sum(widths)-widths[0]/2. + gap/2.)
    else:
        ax.set_xlim(-1.*widths[0]/2. + gap/2., np.sum(widths)-widths[0]/2. - gap/2.)
    ax.set_ylim(0, yTicks[0][-1])#np.max(data_stack))

    # labels
    if xlabel != '':
        plt.xlabel(xlabel)
    if ylabel != '':
        plt.ylabel(ylabel)

#########################################
def plot_celltypesignificance(path, plotdf, args):
    import matplotlib.pyplot as plt
    import numpy as np
    import sys
    print ('plotting celltype significance plot')
    plotdf = plotdf[(plotdf['genecluster'] >= 10)]
    if len(plotdf) < 1:
        sys.exit('Not enough celltypes to plot, please look at the sgcelltype.xls')
    if args.subcommand_name == 'exprbased':
        plotdf = plotdf.sort_values('p-val',ascending=True)
    l = np.log2(plotdf['genecluster'].tolist())
    t = plotdf['p-val'].tolist()
    s = range(1, len(plotdf)+1)
    name = plotdf['celltype'].tolist()
    area = [abs((np.log10(x)) * 15) for x in t]
    color = np.random.random(len(t))
    #print area
    plt.scatter(s, l, s=area, c=color, alpha=0.5)
    plt.grid(True, linestyle=':', color='black')

    # draw a thick red hline at y=0 that spans the xrange
    h = plt.axhline(linewidth=1, color='r', y=5, linestyle='--') #y=args.celltypeClusterSize

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


def heatmap_Sigcelltype(df, path, edgecolors='w', log=False):
    '''

    :param df:
    :param edgecolors:
    :param log:
    :return:
    '''
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
    ## seaborn heatmap
    import seaborn as sns
    sns.set()
    plt.figure(figsize=(15, 20))
    g = sns.heatmap(df, vmin=0, annot=True, cmap="YlGnBu")
    plt.savefig(os.path.join(path, 'GCAM_heatmap_SigCelltype_sns.svg'))
    plt.clf()
    plt.close()
