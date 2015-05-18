__author__ = 'peeyush'
from pandas import read_csv

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
    with open(path) as file:
        for gene in file:
            gene = gene.strip()
            geneList.append(gene.lower())
    return geneList

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
    return geneSyn