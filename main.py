import timeit
from PubmedSearch import Fetch_pmids
from PubmedSearch import Occurrence
from PubmedSearch import ReadFiles
from PubmedSearch import SignificanceTesting


__author__ = 'peeyush'
import pandas as pd



#genenames = ['aatf', 'prmt6', 'ski', 'cd44']
genenames = ReadFiles.get_genes('/home/peeyush/Desktop/genes.txt')
primarygene = genenames
subquery = None
synonym = False
organism = 'human'
annDB = ReadFiles.read_database('/home/peeyush/NetBeansProjects/GCAM-1.0/resources')
cellDB = ReadFiles.celltype_DB('/home/peeyush/NetBeansProjects/GCAM-1.0/resources')
cellSyn = ReadFiles.cell_synonym('/home/peeyush/NetBeansProjects/GCAM-1.0/resources')
start = timeit.default_timer()
print 'User provided gene count:', len(genenames)
if synonym:
    geneSyn = ReadFiles.gene_synonym('/home/peeyush/NetBeansProjects/GCAM-1.0/resources', organism)
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
cellOccu.to_csv('/home/peeyush/Desktop/GCAM_python.csv', sep=',', encoding='utf-8', ignore_index=True)
del annDB, cellDB, cellSyn
###### Scale df for heatmap
significanceDF = SignificanceTesting.SignificanceObject(cellOccu)
significanceDF.heatmap_create()
significanceDF.heatmapdf.to_csv('/home/peeyush/Desktop/GCAM_python_final.csv', sep=',', encoding='utf-8', ignore_index=True)
significanceDF.fisher_test()
significanceDF.pvaldf.to_csv('/home/peeyush/Desktop/GCAM_python_final_pval.csv', sep=',', encoding='utf-8', ignore_index=True)
significanceDF.adjpvaldf.to_csv('/home/peeyush/Desktop/GCAM_python_final_adjpval.csv', sep=',', encoding='utf-8', ignore_index=True)
significanceDF.plot_heatmap()
stop = timeit.default_timer()
print 'Total no. of genes: ', len(genenames)
print 'Time elapsed:', stop-start, ' sec'