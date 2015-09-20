__author__ = 'peeyush'
import timeit, os
from GCAM import Occurrence
from GCAM import FilesFolders
from GCAM import SignificanceTesting
from GCAM import ExpressionAnalysis
from GCAM import Previous_genecheck
from GCAM import plots
import pandas as pd

def gcam_analysis(args, resource_path):
    '''
    Main GCAM function.
    :param args:
    :param resource_path:
    :return: cell occurrence dataframe
    '''
    tstart = timeit.default_timer()
    save_location = args.outdir
    key_celltypes = args.key_celltype_list
    genenames = FilesFolders.get_genes(args.path)
    subquery = args.subquery
    synonym = args.synonym
    primarygene = genenames
    organism = args.org
    ### Reading require databases
    print ('Reading required DBs')
    cellSyn = FilesFolders.cell_synonym(resource_path)

    if synonym:
        geneSyn = FilesFolders.gene_synonym(resource_path, organism)
        genenames = Occurrence.gene2synonym(genenames, geneSyn)
        print ('Gene count after synonym:' + str(len(genenames)))
    occuDF = Previous_genecheck.occurrence_df(genenames, resource_path, subquery)
    cellOccu = Occurrence.joincellsynonym(occuDF, cellSyn)
    if synonym:
        cellOccu = Occurrence.joingenesynonym(cellOccu, primarygene, geneSyn)
    ## Reduced celltypes
    if key_celltypes:
        key_celltypes = FilesFolders.key_celltypes(resource_path)
        cellOccu = cellOccu[cellOccu['celltype'].isin(key_celltypes)]
    #print ('size of new df', len(cellOccu))
    cellOccu = cellOccu.set_index(cellOccu['celltype'])
    cellOccu = cellOccu.drop(['celltype'], axis=1)
    outdir = FilesFolders.create_folders(save_location)
    pd.DataFrame(genenames, columns=['GeneNames']).to_csv(outdir + os.path.sep + 'input_gene_list.csv', sep=',', encoding='utf-8', index=False)
    cellOccu.to_csv(outdir + os.path.sep + 'GCAM_python_occurrence.csv', sep='\t', encoding='utf-8', ignore_index=True)
    ###### Scale df for heatmap and do further analysis
    significanceDF = SignificanceTesting.SignificanceObject(cellOccu)
    significanceDF.heatmapdf_create()
    #try:
    significanceDF.plot_heatmap(outdir)
    significanceDF.fisher_occurrence_test()
    write_result(significanceDF, outdir, key_celltypes)
    ###### Expression analysis of celltype
    subcommand = args.subcommand_name
    if subcommand == "exprbased":
        expressiondf = FilesFolders.read_expression_file(args.exppath)
        expObj = ExpressionAnalysis.ExpressionData(expressiondf)
        expObj.celltype_expression(significanceDF.sigCelltypedf, significanceDF.cellgenedf, outdir)

    tstop = timeit.default_timer()
    print ('Total time elapsed: ' + str(tstop - tstart) + ' sec')
    return cellOccu
    #except:
    #    raise Warning("Genes are not significantly enriched for celltypes or number of queries are < 2")


def write_result(significanceDF, outdir, key_celltypes):
    '''
    Print all the output for genebased analysis.
    :param significanceDF:
    :return:
    '''
    #significanceDF.heatmapdf.to_csv(save_location+'/GCAM_output/GCAM_python_final_occurrence.csv', sep=',', encoding='utf-8', ignore_index=True)
    #significanceDF.pvaldf.to_csv(save_location+'/GCAM_output/GCAM_python_final_pval.csv', sep=',', encoding='utf-8', ignore_index=True)
    #significanceDF.adjpvaldf.to_csv(save_location+'/GCAM_output/GCAM_python_final_adjpval.csv', sep=',', encoding='utf-8', ignore_index=True)

    cellgenedf = significanceDF.cellgenedf[significanceDF.cellgenedf['FDR'] < 0.05]
    cellgenedf.sort(['P-val'], ascending=True)
    if len(cellgenedf)>0:cellgenedf.to_csv(outdir + os.path.sep + 'GCAM_sigenes.csv', sep='\t', encoding='utf-8', index=False)
    else: print('No significant genes for celltype')

    sigCelltypedf = significanceDF.sigCelltypedf[significanceDF.sigCelltypedf['FDR'] < 1]
    plots.stack_barplot(sigCelltypedf, outdir, key_celltypes)
    plots.plot_celltypesignificance(outdir, sigCelltypedf)
    sigCelltypedf.sort(['P-val'], ascending=True)
    if len(sigCelltypedf)>0:sigCelltypedf.to_csv(outdir + os.path.sep + 'GCAM_sigCelltypes.csv', sep='\t', encoding='utf-8', index=False)
    else: print('No significant celltypes')