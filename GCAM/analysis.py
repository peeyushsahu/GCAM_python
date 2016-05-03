import timeit, os
from GCAM import Occurrence
from GCAM import FilesFolders
from GCAM import SignificanceTesting
from GCAM import ExpressionAnalysis, ExpressionClustering
from GCAM import Previous_genecheck
from GCAM import plots
import pandas as pd
import logging

__author__ = 'peeyush'


def gcam_analysis(args, resource_path):
    '''
    Main GCAM function.
    :param args:
    :param resource_path:
    :return: cell occurrence dataframe
    '''
    outdir = FilesFolders.create_folders(args.outdir)
    logging.basicConfig(filename=os.path.join(outdir,'GCAM.log'), level=logging.INFO)
    logging.info('Started')
    if not args.dbpath is None:
        if os.path.exists(args.dbpath):
            logging.info('External resource path used.')
            resource_path = args.dbpath
        else:
            logging.error('Resource dir for db does not exist.')
            print('Resource dir for db does not exist, check --dbpath')
    else:
        logging.info('Default resource path used.')
        resource_path = resource_path

    subcommand = args.subcommand_name
    tstart = timeit.default_timer()

    # Reading require databases
    logging.info('Reading required resource dbs..')
    warnings(args)
    write_parameter(args, subcommand, outdir)

    if subcommand == 'genebased':
        genenames = FilesFolders.get_genes(args.path)
        gene_based(args, resource_path, genenames, outdir)

    if subcommand == 'exprbased':
        expressiondf = FilesFolders.read_expression_file(args.exppath)
        pheno_data = FilesFolders.read_pheno_data(args.phenopath)
        ## print pheno_data['phenotype']
        if args.controlsample is not None and args.controlsample not in list(pheno_data['phenotype']):
            logging.error(args.controlsample+': control sample name not found in phenotype list')
            raise KeyError(args.controlsample+": control sample name not found in phenotype list.")
        genenames, newexprdf = expr_based(outdir, expressiondf, pheno_data, args)
        #print 'genes', genenames
        significance_Df = gene_based(args, resource_path, genenames, outdir)
        plotdf = ExpressionClustering.exprdf4plot(significance_Df, newexprdf, pheno_data, args,
                                                  control=args.controlsample, path=outdir, clusterSize=int(args.celltypeClusterSize)) #'WT_LIV_Abdulah'
        #plots.stack_barplot(plotdf, outdir, key_celltypes=args.key_celltype_list, method=subcommand)
    tstop = timeit.default_timer()
    print ('Total time elapsed: ' + str(tstop - tstart) + ' sec')
    logging.info('Total time elapsed: ' + str(tstop - tstart) + ' sec')
    logging.info('Finished')


def write_parameter(args, subcommand, outdir):
    '''
    Write parameters in a txt file.
    :param args:
    :param subcommand:
    :param outdir:
    :return:
    '''
    if subcommand == 'exprbased':
        parameter = open(os.path.join(outdir, 'parameter.txt'), 'w')
        parameter.write("Analysis type: "+subcommand+"\n")
        parameter.write("Som grid size: "+str(args.som_gridsize)+"\n")
        parameter.write("Minimum no of genes for cell fraction analysis: "+str(args.celltypeClusterSize)+"\n")
        parameter.write("Regression method used: nuSVR \n")
        parameter.write("Consider mean of all samples as reference: "+str(args.meanAsControl)+"\n")
        if not args.meanAsControl:
            parameter.write("Control sample name: "+args.controlsample+"\n")
        parameter.close()


def warnings(args):
    '''
    Gives warning based on paramenter collection.
    :return:
    '''
    subcommand=args.subcommand_name
    if subcommand == 'exprbased':
        if args.controlsample is None:
            if not args.meanAsControl:
                logging.error("Control sample name is not selected.")
                raise ValueError("Please specify control sample name OR set --meanAsControl, -m")
        userCelltype = args.selectCelltypes
        if userCelltype is not None:
            userCelltype = [x.strip() for x in userCelltype.split(',')]
            if len(userCelltype) > 26:
                raise ValueError("Selected cell-types for analysis should be less than 26 and more than 0.")


def gene_based(args, resource_path, genenames, outdir):
    synonym = args.synonym
    organism = args.org
    genenames = genenames
    subquery = args.subquery
    primarygene = genenames
    cellSyn = FilesFolders.cell_synonym(resource_path)
    binom_prob = FilesFolders.read_binom_prob(resource_path)
    pd.DataFrame(genenames, columns=['genenames']).to_csv(outdir + os.path.sep + 'input_gene_list.txt', sep='\t', encoding='utf-8', index=False)

    ocstart = timeit.default_timer()
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
    ocstop = timeit.default_timer()
    logging.info("TC in occurrence analysis:"+str(ocstop - ocstart)+'sec')

    # Scale df for heatmap and do further analysis
    significanceDF = SignificanceTesting.SignificanceObject(cellOccu, binom_prob, resource_path, outdir)
    significanceDF.heatmapdf_create()
    significanceDF.plot_heatmap(outdir)
    significanceDF.fisher_occurrence_test()
    significanceDF.data4radarplot()
    write_result(significanceDF, outdir, args)
    return significanceDF


def expr_based(outdir, expressiondf, pheno_data, args):
    # Expression analysis of celltype
    return ExpressionClustering.SOMclustering(expressiondf, pheno_data, outdir, float(args.som_foldifference), iteration=int(args.somiter))


def write_result(significanceDF, outdir, args):
    '''
    Print all the output for genebased analysis.
    :param significanceDF:
    :return:
    '''
    cellgenedf = significanceDF.cellgenedf  # [significanceDF.cellgenedf['p-val'] < 0.05]
    cellgenedf.sort_values('p-val', ascending=True)
    if len(cellgenedf)>0:
        filtereddf = filter_df(cellgenedf)
        filtereddf.to_excel(os.path.join(outdir, 'GCAM_sigenes.xlsx'), index=False)
    else: print('No significant genes for celltype')

    sigCelltypedf = significanceDF.sigCelltypedf[significanceDF.sigCelltypedf['FDR'] < 0.05]
    if len(sigCelltypedf) > 1:
        plots.plot_celltypesignificance(outdir, sigCelltypedf, args)
        sigCelltypedf.sort_values('genecluster', ascending=True)
        sigCelltypedf.to_excel(os.path.join(outdir, 'GCAM_sigCelltypes.xlsx'), index=False)
    else: print('No significant celltypes')


def filter_df(df):
    '''
    Filter sigenes df to print only 3 most significant gene-celltype entries.
    :param df:
    :return:
    '''
    new_df = pd.DataFrame(columns=df.columns)
    df_gr = df.groupby('gene')
    for k, gr in df_gr:
        new_gr = gr.sort_values('p-val', ascending=True)
        rows = 1
        for i, r in new_gr.iterrows():
            if rows <= 3:
                new_df = new_df.append(r)
                rows += 1
    return new_df