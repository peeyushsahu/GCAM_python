__author__ = 'peeyush'
from pandas import read_csv
import os, sys, stat
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
    print ('Output directory created: ' + os.path.join(path, 'GCAM_output_'+str(time.strftime("%d_%m_%Y"))+'_'+str(time.strftime("%H_%M_%S"))))
    npath = os.path.join(path, 'GCAM_output_'+str(time.strftime("%d_%m_%Y"))+'_'+str(time.strftime("%H_%M_%S")))
    if not os.path.exists(npath):
        os.makedirs(npath)
        #print('Changed permission')
        os.chmod(npath, 0o777)
        #os.chmod(npath, mode=777)
    return npath

def read_annotation_database(path):
    '''
    Reading annotation DB as pandas dataframe
    :param path:
    :return:
    '''
    try:
        annoDB = read_csv(path + os.path.sep + 'pmid_celltype_index_final.txt', header=0, sep="\t")
    except:
        raise ValueError("annotation db doesnot exist.")
    annoDB = annoDB.set_index(['pmid'])
    return annoDB


def celltype_DB(path):
    '''
    Import celltype database.
    :param path:
    :return:
    '''
    try:
        cellDB = read_csv(path + os.path.sep + 'cellTypes.csv', header=None, sep=',')
    except:
        raise ValueError("celltype db does not exist.")
    cellDB.columns = ['celltype']
    cellDB['celltype'] = cellDB['celltype'].str.lower()
    #print cellDB
    return cellDB


def cell2tissue_DB(path):
    '''
    Import celltype database.
    :param path:
    :return:
    '''
    try:
        cell2tissue = read_csv(os.path.join(path, 'cell2tissue.txt'), header=0, sep='\t', index_col=0)
    except:
        raise ValueError("cell2tissue db does not exist.")
    #print cellDB
    return cell2tissue


def cell_synonym(path):
    '''
    Import cell synonym database.
    :param path:
    :return:
    '''
    try:
        cellSyn = read_csv(path + os.path.sep + 'cell_type_synonyms_python.csv', header=0, sep=',')
    except:
        raise ValueError("Synonym db does not exist.")
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
                if len(gene) > 0:
                    geneList.append(gene.lower())
    except IOError:
        print ("Error: Query file does not appear to exist.")
        sys.exit(1)
    f_genelist = list(set(geneList))
    print ('Size of user provided gene list:' + str(len(geneList)))
    print ('No. of genes after removing duplicates:' + str(len(f_genelist)))
    if len(f_genelist) < 5:
        raise IOError("Minimum gene requirement for the analysis is 5.")
    return f_genelist

def key_celltypes(path):
    '''
    Read input gene list.
    :param path:
    :return:
    '''
    celltypeList = []
    path = path + os.path.sep + 'key_cellTypes.csv'
    try:
        with open(path) as file:
            for cells in file:
                cells = cells.strip()
                celltypeList.append(cells.lower())
    except IOError:
        print ("Error: Key celltype db does not exist.")
        sys.exit(1)
    f_celltypes = list(set(celltypeList))
    print ('Size of key gene list:' + str(len(celltypeList)))
    #print 'No. of genes after removing duplicates:', len(f_genelist)
    return f_celltypes

def gene_synonym(path, organism):
    '''
    Reads synonym file for genes
    :param path:
    :param organism:
    :return:
    '''
    try:
        if organism == 'human':
            geneSyn = read_csv(path + os.path.sep + 'Human_synonym.txt', header=None, sep='\t')
        elif organism == 'mouse':
            geneSyn = read_csv(path + os.path.sep + 'Mouse_synonym.txt', header=None, sep='\t')
    except:
        raise ValueError("Synonym db does not exist.")
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
        # print ('exprs path:'+ path)
        expressiondf = read_csv(path, header=0, sep="\t", index_col=0)
        # print(expressiondf.head(1))
        expressiondf.index = expressiondf.index.to_series().str.lower()
        expressiondf = expressiondf[expressiondf.index.to_series().str.len() >= 3] ##Removing gene names with only 2 letters.
    except IOError:
        print ("Error: Expression File does not appear to exist.")
        sys.exit(0)
    return expressiondf


def read_previous_occurrence_table(resource_path):
    '''
    This will read the occurrence database for already analysed genes to save time.
    :param resource_path:
    :return:
    '''
    try:
        print ('reading previously analysed genes')
        gene_occu_db = read_csv(os.path.join(resource_path, 'gene_occu_db.csv'), header=0, sep=",", index_col=0)
    except:
        print ("Warning: Creating Gene occurrence db, analysis will take longer.")
        return None, False
    return gene_occu_db, True


def read_binom_prob(resource_path):
    '''
    This will read the occurrence database for already analysed genes to save time.
    :param resource_path:
    :return:
    '''
    try:
        #print ('binomial probability calculation...')
        cell_binom = read_csv(os.path.join(resource_path, 'cellType_binom_prob.txt'), header=0, sep="\t")
    except:
        raise ValueError('Binomial probabilities not found')
    return cell_binom


def read_pheno_data(path):
    '''
    Reads pheno-data for analysis
    :param path:
    :return:
    '''
    try:
        #print ('pheno path:'+ path)
        phenodf = read_csv(path, header=0, sep="\t")
        if not 'phenotype' in phenodf.columns or not 'sample' in phenodf.columns:
            print ('Error phenotype: please name columns as phenotype and sample')
            sys.exit(0)
    except IOError:
        print ("Error: pheno File does not appear to exist.")
        sys.exit(0)
    return phenodf

def read_cell_subtractdf(path):
    '''
    Reads pheno-data for analysis
    :param path:
    :return:
    '''
    try:
        cellsub_df = read_csv(os.path.join(path, 'celltype_subtract.txt'), header=0, sep="\t")
        #print cellsub_df
    except IOError:
        print ("Error: cell_subtract db does not appear to exist.")
        sys.exit(0)
    return cellsub_df