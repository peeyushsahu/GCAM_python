import timeit
from PubmedSearch import Fetch_pmids
from PubmedSearch import Occurrence
from PubmedSearch import FilesFolders
from PubmedSearch import SignificanceTesting
from PubmedSearch import ExpressionAnalysis


__author__ = 'peeyush'



start = timeit.default_timer()
save_location = '/home/peeyush/Desktop'
FilesFolders.create_folders(save_location)
genenames = ['aatf', 'prmt6', 'ski', 'cd44']
#genenames = FilesFolders.get_genes('/home/peeyush/Desktop/genes.txt')
primarygene = genenames
subquery = None
synonym = False
organism = 'human'
annDB = FilesFolders.read_database('/home/peeyush/NetBeansProjects/GCAM-1.0/resources')
cellDB = FilesFolders.celltype_DB('/home/peeyush/NetBeansProjects/GCAM-1.0/resources')
cellSyn = FilesFolders.cell_synonym('/home/peeyush/NetBeansProjects/GCAM-1.0/resources')
expressiondf = FilesFolders.read_expression_file('/home/peeyush/Desktop/GCAM_expression.csv')


if synonym:
    geneSyn = FilesFolders.gene_synonym('/home/peeyush/NetBeansProjects/GCAM-1.0/resources', organism)
    genenames = Occurrence.gene2synonym(genenames, geneSyn)
    print 'Gene count after synonym:', len(genenames)

occuDF = cellDB
for i in genenames:
    gene = Fetch_pmids.Genes(i)
    gene.get_pmids()
    gene.get_pmid_pos(annoDB=annDB)
    occuDF = gene.get_occurrence(cellDB=occuDF)

cellOccu = Occurrence.joincellsynonym(occuDF, cellSyn)
if synonym:
    cellOccu = Occurrence.joingenesynonym(cellOccu, primarygene, geneSyn)
cellOccu = cellOccu.set_index(cellOccu['celltype'])
cellOccu = cellOccu.drop(['celltype'], axis=1)
cellOccu.to_csv(save_location+'/GCAM_output/GCAM_python_occurrence.csv', sep=',', encoding='utf-8', ignore_index=True)
del annDB, cellDB, cellSyn

###### Scale df for heatmap
significanceDF = SignificanceTesting.SignificanceObject(cellOccu)
significanceDF.heatmapdf_create()
significanceDF.heatmapdf.to_csv(save_location+'/GCAM_output/GCAM_python_final_occurrence.csv', sep=',', encoding='utf-8', ignore_index=True)
significanceDF.fisher_test()
significanceDF.pvaldf.to_csv(save_location+'/GCAM_output/GCAM_python_final_pval.csv', sep=',', encoding='utf-8', ignore_index=True)
significanceDF.adjpvaldf.to_csv(save_location+'/GCAM_output/GCAM_python_final_adjpval.csv', sep=',', encoding='utf-8', ignore_index=True)
significanceDF.cellgenedf.to_csv(save_location+'/GCAM_output/GCAM_python_final_cellGene.csv', sep=',', encoding='utf-8', ignore_index=True)
significanceDF.sigCelltypedf.to_csv(save_location+'/GCAM_output/GCAM_python_final_SigCelltypes.csv', sep=',', encoding='utf-8', ignore_index=True)
significanceDF.plot_heatmap(save_location)


###### Expression analysis of celltype
expObj = ExpressionAnalysis.ExpressionData(expressiondf)
expObj.celltype_expression(significanceDF.sigCelltypedf, significanceDF.cellgenedf, save_location)
expObj.plotdf.to_csv(save_location+'/GCAM_output/GCAM_python_final_celltype_vs_expression.csv', sep=',', encoding='utf-8', ignore_index=True)

stop = timeit.default_timer()
print 'Total no. of genes: ', len(genenames)
print 'Time elapsed:', stop-start, ' sec'
