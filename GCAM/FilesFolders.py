__author__ = 'peeyush'
from pandas import read_csv
import os, sys

def create_folders(path):
    '''
    folders = ["overlap",
               "differential",
               "filtered",
               "plots",
               "seq4motif",
               "valley_peaks",
               "CpG",
               "density_based_motif"]
    for folder in folders:
    '''
    print 'Output directory created: ' + path +'/GCAM_output'
    npath = path +'/GCAM_output'
    if not os.path.exists(npath):
        os.makedirs(npath)

def read_database(path):
    '''
    Reading annotation DB as pandas dataframe
    :param path:
    :return:
    '''
    annoDB = read_csv(path + '/pmid_celltype_index_final.txt', header=0, sep="\t")
    annoDB = annoDB.set_index(['pmid'])
    return annoDB


def celltype_DB(path):
    '''
    Import celltype database.
    :param path:
    :return:
    '''
    cellDB = read_csv(path + '/cellTypes.csv', header=None, sep=',')
    cellDB.columns = ['celltype']
    return cellDB


def cell_synonym(path):
    '''
    Import cell synonym database.
    :param path:
    :return:
    '''
    cellSyn = read_csv(path + '/cell_type_synonyms_python.csv', header=0, sep=',')
    return cellSyn

def get_genes(path):
    '''
    Read input gene list.
    :param path:
    :return:
    '''
    geneList = []
    try:
        with open(path) as file:
            for gene in file:
                gene = gene.strip()
                geneList.append(gene.lower())
    except IOError:
        print "Error: File does not appear to exist."
        sys.exit(1)
    f_genelist = list(set(geneList))
    print 'Size of user provided gene list:', len(geneList)
    print 'No. of genes after remove duplicates:', len(f_genelist)
    return f_genelist

def gene_synonym(path, organism):
    '''
    Reads synonym file for genes
    :param path:
    :param organism:
    :return:
    '''
    if organism == 'human':
        geneSyn = read_csv(path + '/Human_synonym.txt', header=None, sep='\t')
    elif organism == 'mouse':
        geneSyn = read_csv(path + '/Mouse_synonym.txt', header=None, sep='\t')
    geneSyn.columns = ['gene', 'synonym']
    geneSyn['gene'] = geneSyn['gene'].str.lower()
    geneSyn = geneSyn.set_index(geneSyn['gene'])
    return geneSyn

def read_expression_file(path):
    '''
    Reads expression data for analysis
    :param path:
    :return:
    '''
    expressiondf = read_csv(path, header=0, sep=",")
    expressiondf['genes'] = expressiondf['genes'].str.lower()
    expressiondf = expressiondf.set_index(expressiondf['genes'])
    return expressiondf