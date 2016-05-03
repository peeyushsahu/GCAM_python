__author__ = 'peeyush'
import pandas as pd
import os


class ExpressionData():

    def __init__(self, expdf):
        self.expressiondf = expdf
        self.plotdf = None

    def celltype_expression(self, sigCelltypedf, gene2celltypedf, path):
        gene2cell_group = gene2celltypedf.groupby(gene2celltypedf['CellType'])
        expression_df = self.expressiondf
        #print expression_df
        fold_change_column = []
        for i in expression_df.columns.tolist():
            if 'Fold' in i:fold_change_column.append(i)
        for column in fold_change_column:
            plotdf = pd.DataFrame()
            for k, v in sigCelltypedf.iterrows():
                if v['P-val'] < 0.1:
                    #print v['CellType']
                    expression_list = []
                    df = gene2cell_group.get_group(v['celltype'])
                    #print df.shape
                    for key, val in df.iterrows():
                        if val['P-val'] < 0.001:
                            expression = expression_df.loc[val['Genes'], column]
                            #print 'expression for '+val['Genes'], expression
                            if type(expression) == pd.Series:expression_list.append(max(expression))
                            else: expression_list.append(expression)
                    #print 'exprssion list for gene', v['CellType'], expression_list
                    fold_change = 0
                    if len(expression_list) > 0:fold_change = sum(expression_list)/len(expression_list)
                    plotdf = plotdf.append(pd.Series([v['celltype'], v['genecluster'],v['P-val'], fold_change]), ignore_index=True)
                    #print [v['CellType'], v['P-val'], fold_change]
            if len(plotdf) > 0:
                plotdf.columns = ['celltype', 'genecluster', 'p-val', 'relative_expression']
                plotdf = plotdf[plotdf['p-val'] < 0.5]
                plotdf.to_csv(path + os.path.sep + 'GCAM_python_final_celltype_vs_expression_'+column+'.csv', sep=',', encoding='utf-8', ignore_index=True)
                plot_expressionvseignificance(path, plotdf, column)
            else:
                print ('No significant celltypes for the data')


def plot_expressionvseignificance(path, plotdf, column):
    import matplotlib.pyplot as plt
    import numpy as np
    import math
    plotdf = plotdf[plotdf['p-val'] < 0.05]
    print ('plotting significance plot')
    plotdf = plotdf.sort(['p-val'], ascending=True)
    l = plotdf['genecluster'].tolist()
    t = plotdf['p-val'].tolist()
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
    plt.ylim(min(s)-2,max(s)+2)
    plt.title('Relative Cell-type expression plot' + column, fontsize=14)
    plt.xlabel('Gene cluster size', fontsize=12)
    plt.ylabel('Average Fold Change', fontsize=12)
    #plt.tight_layout()
    plt.savefig(path + os.path.sep + column + 'GCAM_celltype_VS_expresiion.png')
    plt.clf()

