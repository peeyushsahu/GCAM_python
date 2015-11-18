__author__ = 'peeyush'
import timeit, os
from GCAM import Occurrence
from GCAM import FilesFolders
from GCAM import SignificanceTesting
from GCAM import ExpressionAnalysis, ExpressionClustering
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
    subcommand = args.subcommand_name
    tstart = timeit.default_timer()
    ### Reading require databases
    print ('Reading required DBs')
    outdir = FilesFolders.create_folders(args.outdir)
    if subcommand == 'exprbased':
        expressiondf = FilesFolders.read_expression_file(args.exppath)
        pheno_data = FilesFolders.read_pheno_data(args.phenopath)
    #print 'Keycelltype', args.key_celltype_list, args.outdir, args.synonym
    if subcommand == 'genebased':
        genenames = FilesFolders.get_genes(args.path)
        gene_based(args, resource_path, genenames, outdir)
    if subcommand == 'exprbased':
        genenames, newexprdf = expr_based(outdir, expressiondf, pheno_data, args)
        #print 'genes', genenames
        significance_Df = gene_based(args, resource_path, genenames, outdir)
        plotdf = ExpressionClustering.exprdf4plot(significance_Df, newexprdf, pheno_data, args.regMethod,
                                                  control=args.controlsample, path=outdir, clusterSize=int(args.celltypeClusterSize)) #'WT_LIV_Abdulah'
        plots.stack_barplot(plotdf, outdir, key_celltypes=args.key_celltype_list, method=subcommand)
    tstop = timeit.default_timer()
    print ('Total time elapsed: ' + str(tstop - tstart) + ' sec')


def gene_based(args, resource_path, genenames, outdir):
    synonym = args.synonym
    organism = args.org
    genenames = genenames
    subquery = args.subquery
    primarygene = genenames
    cellSyn = FilesFolders.cell_synonym(resource_path)
    pd.DataFrame(genenames, columns=['GeneNames']).to_csv(outdir + os.path.sep + 'input_gene_list.csv', sep=',', encoding='utf-8', index=False)
    if synonym:
        geneSyn = FilesFolders.gene_synonym(resource_path, organism)
        genenames = Occurrence.gene2synonym(genenames, geneSyn)
        print ('Gene count after synonym:' + str(len(genenames)))
    occuDF = Previous_genecheck.occurrence_df(genenames, resource_path, subquery)
    cellOccu = Occurrence.joincellsynonym(occuDF, cellSyn)
    if synonym:
        cellOccu = Occurrence.joingenesynonym(cellOccu, primarygene, geneSyn)
    # Reduced celltypes
    if args.key_celltype_list:
        key_celltypes = FilesFolders.key_celltypes(resource_path)
        cellOccu = cellOccu[cellOccu['celltype'].isin(key_celltypes)]
    # print ('size of new df', len(cellOccu))
    cellOccu = cellOccu.set_index(cellOccu['celltype'])
    cellOccu = cellOccu.drop(['celltype'], axis=1)
    ## Subtract eg. cd4 t cell from t cell
    #cellOccu = Occurrence.subtract_cellnamepeat(cellOccu, resource_path)
    cellOccu.to_csv(outdir + os.path.sep + 'GCAM_python_occurrence.csv', sep='\t', encoding='utf-8', ignore_index=True)
    # Scale df for heatmap and do further analysis
    significanceDF = SignificanceTesting.SignificanceObject(cellOccu)
    significanceDF.heatmapdf_create()
    significanceDF.plot_heatmap(outdir)
    significanceDF.fisher_occurrence_test()
    write_result(significanceDF, outdir, key_celltypes=args.key_celltype_list)
    return significanceDF


def expr_based(outdir, expressiondf, pheno_data, args):
    # Expression analysis of celltype
    return ExpressionClustering.SOMclustering(expressiondf, pheno_data, outdir, float(args.som_foldifference), iteration=int(args.somiter))

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