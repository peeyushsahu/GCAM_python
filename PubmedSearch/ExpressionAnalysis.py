__author__ = 'peeyush'
import numpy as np
import pandas as pd


class ExpressionData():

    def __init__(self, expdf):
        self.expressiondf = expdf
        self.plotdf = None

    def celltype_expression(self, sigCelltypedf, gene2celltypedf, path):
        plotdf = pd.DataFrame()
        gene2cell_group = gene2celltypedf.groupby(gene2celltypedf['CellType'])
        for k, v in sigCelltypedf.iterrows():
            if v['P-val'] < 1:
                #print v['CellType']
                expression_list = []
                df = gene2cell_group.get_group(v['CellType'])
                #print df.shape
                for key, val in df.iterrows():
                    expression = self.expressiondf.loc[val['Genes'], 'FoldChange']
                    expression_list.append(expression)
                plotdf = plotdf.append(pd.Series([v['CellType'], v['P-val'], sum(expression_list)/len(expression_list)]), ignore_index=True)
                #print [v['CellType'], v['P-val'], sum(expression_list)/len(expression_list)]
        #print plotdf
        plotdf.columns = ['celltype', 'p-val', 'expression']
        self.plotdf = plotdf
        ## Function for plot
        self.plot_expressionvseignificance(path)


    def plot_expressionvseignificance(self, path):
        import matplotlib.pyplot as plt
        print 'plotting significance plot'
        plotdf = self.plotdf
        print plotdf

        t = plotdf['p-val'].tolist()
        s = plotdf['expression'].tolist()
        name = plotdf['celltype'].tolist()
        #t = np.arange(0.0000001, 0.2, 0.005)
        #s = 10 * np.random.random_sample((len(t),))-5

        plt.plot(t, s, 'o', markersize=5, markerfacecolor='black')
        # draw a thick red hline at y=0 that spans the xrange
        l = plt.axhline(linewidth=1, color='black')

        # draw a default vline at x=1 that spans the yrange
        l = plt.axvline(linewidth=1, color='black', x=0.05)

        # draw color background parition on y axis
        p = plt.axhspan(0, 10, facecolor='r', alpha=0.5)
        p = plt.axhspan(0, -10, facecolor='b', alpha=0.5)

        for i in range(0,len(t)):
            if t[i] < 0.05:
                plt.annotate(name[i], xy=(t[i], s[i]),  xycoords='data',
                            xytext=(30, 30), textcoords='offset points',
                            arrowprops=dict(arrowstyle="->",
                                            connectionstyle="arc3,rad=0"))
        plt.xlim(min(t)-0.01, max(t)+0.01)
        plt.title('Expression vs Significance')
        plt.xlabel('P-Value')
        plt.ylabel('Average Fold Change')
        #plt.show()
        plt.savefig(path+'/GCAM_output/GCAM_celltype_VS_expresiion.png')
        plt.close()