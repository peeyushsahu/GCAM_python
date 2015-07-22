__author__ = 'peeyush'
import timeit, os
from GCAM import Occurrence
from GCAM import FilesFolders
from GCAM import SignificanceTesting
from GCAM import ExpressionAnalysis
from GCAM import Previous_genecheck

def gcam_analysis(args, resource_path):
    '''
    Main GCAM function.
    :param args:
    :param resource_path:
    :return: cell occurrence dataframe
    '''
    tstart = timeit.default_timer()
    save_location = args.outdir
    key_celltypes = True
    genenames = FilesFolders.get_genes(args.path)
    subquery = args.subquery
    synonym = args.synonym
    primarygene = genenames
    organism = args.org
    ### Reading require databases
    print 'Reading required DBs'
    cellSyn = FilesFolders.cell_synonym(resource_path)

    if synonym:
        geneSyn = FilesFolders.gene_synonym(resource_path, organism)
        genenames = Occurrence.gene2synonym(genenames, geneSyn)
        print 'Gene count after synonym:', len(genenames)
    occuDF = Previous_genecheck.occurrence_df(genenames, resource_path, subquery)
    cellOccu = Occurrence.joincellsynonym(occuDF, cellSyn)
    if synonym:
        cellOccu = Occurrence.joingenesynonym(cellOccu, primarygene, geneSyn)
    ## Reduced celltypes
    if key_celltypes:
        key_celltypes = FilesFolders.key_celltypes(resource_path)
        cellOccu = cellOccu[cellOccu['celltype'].isin(key_celltypes)]
    print 'size of new df', len(cellOccu)
    cellOccu = cellOccu.set_index(cellOccu['celltype'])
    cellOccu = cellOccu.drop(['celltype'], axis=1)
    outdir = FilesFolders.create_folders(save_location)
    cellOccu.to_csv(outdir + os.path.sep + 'GCAM_python_occurrence.csv', sep=',', encoding='utf-8', ignore_index=True)
    ###### Scale df for heatmap and do further analysis
    significanceDF = SignificanceTesting.SignificanceObject(cellOccu)
    significanceDF.heatmapdf_create()
    significanceDF.fisher_occurrence_test()
    significanceDF.plot_heatmap(outdir)
    write_result(significanceDF, outdir)
    ###### Expression analysis of celltype
    subcommand = args.subcommand_name
    if subcommand == "exprbased":
        expressiondf = FilesFolders.read_expression_file(args.exppath)
        expObj = ExpressionAnalysis.ExpressionData(expressiondf)
        expObj.celltype_expression(significanceDF.sigCelltypedf, significanceDF.cellgenedf, outdir)

    tstop = timeit.default_timer()
    print 'Total time elapsed: ', (tstop - tstart), ' sec'
    return cellOccu


def write_result(significanceDF, outdir):
    '''
    Print all the output for genebased analysis.
    :param significanceDF:
    :return:
    '''
    #significanceDF.heatmapdf.to_csv(save_location+'/GCAM_output/GCAM_python_final_occurrence.csv', sep=',', encoding='utf-8', ignore_index=True)
    #significanceDF.pvaldf.to_csv(save_location+'/GCAM_output/GCAM_python_final_pval.csv', sep=',', encoding='utf-8', ignore_index=True)
    #significanceDF.adjpvaldf.to_csv(save_location+'/GCAM_output/GCAM_python_final_adjpval.csv', sep=',', encoding='utf-8', ignore_index=True)
    cellgenedf = significanceDF.cellgenedf[significanceDF.cellgenedf['P-val'] < 0.05]
    cellgenedf.to_csv(outdir + os.path.sep + 'GCAM_python_final_cellGene.csv', sep=',', encoding='utf-8', ignore_index=True)
    sigCelltypedf = significanceDF.sigCelltypedf[significanceDF.sigCelltypedf['P-val'] < 0.05]
    sigCelltypedf.to_csv(outdir + os.path.sep + 'GCAM_python_final_SigCelltypes.csv', sep=',', encoding='utf-8', ignore_index=True)