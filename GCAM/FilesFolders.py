__author__ = 'peeyush'
from pandas import read_csv
import os, sys
import time

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
    print 'Output directory created: ' + path +'/GCAM_output_'+str(time.strftime("%d/%m/%Y"))+str(time.strftime("%H:%M:%S"))
    npath = path + os.path.sep +'GCAM_output_'+str(time.strftime("%d-%m-%Y"))+'_'+str(time.strftime("%H:%M:%S"))
    if not os.path.exists(npath):
        os.makedirs(npath)
        #os.chmod(npath, mode=777)
    return npath

def read_database(path):
    '''
    Reading annotation DB as pandas dataframe
    :param path:
    :return:
    '''
    annoDB = read_csv(path + os.path.sep + 'pmid_celltype_index_final.txt', header=0, sep="\t")
    annoDB = annoDB.set_index(['pmid'])
    return annoDB


def celltype_DB(path):
    '''
    Import celltype database.
    :param path:
    :return:
    '''
    cellDB = read_csv(path + os.path.sep + 'cellTypes.csv', header=None, sep=',')
    cellDB.columns = ['celltype']
    cellDB['celltype'] = cellDB['celltype'].str.lower()
    return cellDB


def cell_synonym(path):
    '''
    Import cell synonym database.
    :param path:
    :return:
    '''
    cellSyn = read_csv(path + os.path.sep + 'cell_type_synonyms_python.csv', header=0, sep=',')
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
    print 'No. of genes after removing duplicates:', len(f_genelist)
    return f_genelist

def gene_synonym(path, organism):
    '''
    Reads synonym file for genes
    :param path:
    :param organism:
    :return:
    '''
    if organism == 'human':
        geneSyn = read_csv(path + os.path.sep + 'Human_synonym.txt', header=None, sep='\t')
    elif organism == 'mouse':
        geneSyn = read_csv(path + os.path.sep + 'Mouse_synonym.txt', header=None, sep='\t')
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
    try:
        print 'exprs path:',path
        expressiondf = read_csv(path, header=0, sep=",")
        if 'SYMBOL' in expressiondf.columns:
            expressiondf['SYMBOL'] = expressiondf['SYMBOL'].str.lower()
            expressiondf = expressiondf.set_index(expressiondf['SYMBOL'])
        else:
            print 'Error: please name columns as SYMBOL and FoldChange'
            sys.exit(0)
    except IOError:
        print "Error: Expression File does not appear to exist."
        sys.exit(0)
    return expressiondf

def read_previous_occurrence_table(resource_path):
    '''
    This will read the occurrence database for already analysed genes to save time.
    :param resource_path:
    :return:
    '''
    try:
        print 'reading previously analysed genes'
        gene_occu_db = read_csv(resource_path + os.path.sep + 'gene_occu_db.csv', header=0, sep=",", index_col=0)
    except:
        print "Warning: gene_occu_db does not appear to exist. Analysis could take more time."
        return None, False
    return gene_occu_db, True