__author__ = 'peeyush'
import timeit, os
from GCAM import Fetch_pmids
from GCAM import Occurrence
from GCAM import FilesFolders
from GCAM import SignificanceTesting
from GCAM import ExpressionAnalysis

def gcam_analysis(args):

    tstart = timeit.default_timer()

    #save_location = '/home/peeyush/Desktop'
    save_location = args.outdir
    FilesFolders.create_folders(save_location)

    #genenames = ['aatf', 'prmt6', 'ski']
    genenames = FilesFolders.get_genes(args.path)

    subquery = args.subquery
    synonym = args.synonym

    genenames = genenames[:50]
    primarygene = genenames
    organism = args.org

    database_path = databasepath()
    database_path = os.path.sep.join(database_path.strip().split(os.path.sep)[:-2])
    database_path = os.path.sep.join([database_path, "GCAM", "resources"])
    print database_path
    annDB = FilesFolders.read_database(database_path)   #'/home/peeyush/NetBeansProjects/GCAM-1.0/resources'
    cellDB = FilesFolders.celltype_DB(database_path)
    cellSyn = FilesFolders.cell_synonym(database_path)

    if synonym:
        geneSyn = FilesFolders.gene_synonym(database_path, organism)
        genenames = Occurrence.gene2synonym(genenames, geneSyn)
        print 'Gene count after synonym:', len(genenames)
    occuDF = cellDB
    fetch_time = 0
    occu_time = 0
    total_abstract = 0
    abs_in_DB = 0
    for gene in genenames:
        #print gene
        fstart = timeit.default_timer()
        GeneObj = Fetch_pmids.Genes(gene)
        GeneObj.get_pmids()
        fstop = timeit.default_timer()
        fetch_time = fetch_time + (fstop - fstart)
        total_abstract += len(GeneObj.pmids) # calcuate total no of abstracts

        ostart = timeit.default_timer()
        GeneObj.get_pmid_pos(annoDB=annDB)
        abs_in_DB += len(GeneObj.cellinpmid)
        occuDF = GeneObj.get_occurrence(cellDB=occuDF)
        ostop = timeit.default_timer()
        occu_time = occu_time + (ostop - ostart)

    cellOccu = Occurrence.joincellsynonym(occuDF, cellSyn)
    if synonym:
        cellOccu = Occurrence.joingenesynonym(cellOccu, primarygene, geneSyn)
    cellOccu = cellOccu.set_index(cellOccu['celltype'])
    cellOccu = cellOccu.drop(['celltype'], axis=1)
    cellOccu.to_csv(save_location+'/GCAM_output/GCAM_python_occurrence.csv', sep=',', encoding='utf-8', ignore_index=True)
    #del annDB, cellDB, cellSyn

    ###### Scale df for heatmap
    significanceDF = SignificanceTesting.SignificanceObject(cellOccu)
    significanceDF.heatmapdf_create()
    significanceDF.heatmapdf.to_csv(save_location+'/GCAM_output/GCAM_python_final_occurrence.csv', sep=',', encoding='utf-8', ignore_index=True)
    significanceDF.fisher_occurrence_test()
    significanceDF.pvaldf.to_csv(save_location+'/GCAM_output/GCAM_python_final_pval.csv', sep=',', encoding='utf-8', ignore_index=True)
    significanceDF.adjpvaldf.to_csv(save_location+'/GCAM_output/GCAM_python_final_adjpval.csv', sep=',', encoding='utf-8', ignore_index=True)
    significanceDF.cellgenedf.to_csv(save_location+'/GCAM_output/GCAM_python_final_cellGene.csv', sep=',', encoding='utf-8', ignore_index=True)
    significanceDF.sigCelltypedf.to_csv(save_location+'/GCAM_output/GCAM_python_final_SigCelltypes.csv', sep=',', encoding='utf-8', ignore_index=True)
    significanceDF.plot_heatmap(save_location)


    ###### Expression analysis of celltype
    subcommand = args.subcommand_name
    if subcommand == "exprbased":
        expressiondf = FilesFolders.read_expression_file('/home/peeyush/Desktop/GCAM_expression.csv')

        expObj = ExpressionAnalysis.ExpressionData(expressiondf)
        expObj.celltype_expression(significanceDF.sigCelltypedf, significanceDF.cellgenedf, save_location)
        expObj.plotdf.to_csv(save_location+'/GCAM_output/GCAM_python_final_celltype_vs_expression.csv', sep=',', encoding='utf-8', ignore_index=True)

    tstop = timeit.default_timer()
    print 'Total no. of genes: ', len(genenames)
    print 'Total no. of abstarcts: ', total_abstract
    print 'Total no. of abstarcts annotated in DB:', abs_in_DB
    print 'Total time elapsed: ', (tstop - tstart), ' sec'
    print 'Total time for pmid fetch: ', fetch_time, ' sec'
    print 'Total time for occurrence analysis: ', occu_time, ' sec'


def databasepath():
    try:
        root = __file__
        if os.path.islink(root):
            root = os.path.realpath(root)
            print os.path.dirname(os.path.abspath(root))
        return os.path.dirname(os.path.abspath(root))
    except:
        print "I'm sorry, but something is wrong."
        print "There is no __file__ variable. Please contact the author."
        sys.exit()