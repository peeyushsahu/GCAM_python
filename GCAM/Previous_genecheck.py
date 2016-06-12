__author__ = 'peeyush'

from GCAM import FilesFolders
import timeit, sys, os
from GCAM import Fetch_pmids
import pandas as pd


def check_database_update(annoDB, cellDB, resource_path):
    '''
    This function checks for update in annotation db and celltype db.
    If found then delete gene_occu_db.
    :param annoDB:
    :param cellDB:
    :return:
    '''
    annoDB_size = len(annoDB)
    cellDB_size = len(cellDB)
    size = []
    try:
        file = open(resource_path + os.path.sep + 'db_version', 'r')
        for line in file:
            size.append(int(line.split(':')[1]))
            #print 'check_database', int(line.split(':')[1])
        file.close()
        if annoDB_size not in size or cellDB_size not in size:
            os.remove(resource_path + os.path.sep + 'gene_occu_db.csv')
            file = open(resource_path + os.path.sep + 'db_version', 'w')
            lines = ['annoDB:' + str(annoDB_size), '\ncellDB:' + str(cellDB_size)]
            file.writelines(lines)
            file.close()
    except:
        file = open(resource_path + os.path.sep + 'db_version', 'w')
        lines = ['annoDB:' + str(annoDB_size), '\ncellDB:' + str(cellDB_size)]
        file.writelines(lines)
        file.close()


def check_old_analysed_genes(genenames, dataframe):
    '''
    This function will check if the gene has already been analysed and retrieve its occurrence.
    :param genenames:
    :param resource_path:
    :return:
    '''
    #dataframe = FilesFolders.read_previous_occurrence_table(resource_path)
    new_genelist = []
    has_genes = []
    if not dataframe is None:
        for gene in genenames:
            if gene not in dataframe.columns:
                new_genelist.append(gene)
            if gene in dataframe.columns:
                has_genes.append(gene)
    foundgenes_df = dataframe[has_genes]
    return new_genelist, foundgenes_df


def occurrence_df(genenames, resource_path, subquery):
    '''
    This function will prepare the occurrence df for genes.
    :param genenames:
    :param resource_path:
    :return:
    '''
    annDB = FilesFolders.read_annotation_database(resource_path)
    cellDB = FilesFolders.celltype_DB(resource_path)
    check_database_update(annDB, cellDB, resource_path)
    print ('Checking for pre analysed genes....')
    dataframe, created = FilesFolders.read_previous_occurrence_table(resource_path)
    join = False
    if dataframe is not None:
        new_genenames, foundgenes_df = check_old_analysed_genes(genenames, dataframe)
        if len(new_genenames) > 0: join = True
    else:
        foundgenes_df = pd.DataFrame()
        new_genenames = genenames
    print ('Reading required DBs')
    occuDF = cellDB
    total_abstract = 0
    abs_in_DB = 0
    count = 0 + len(new_genenames)
    for gene in new_genenames:
        sys.stdout.write("\rGenes remain for analyse:%d" % count)
        sys.stdout.flush()
        #print gene
        GeneObj = Fetch_pmids.Genes(gene=gene, subquery=subquery, resource_path=resource_path)
        GeneObj.get_pmids()
        total_abstract += len(GeneObj.pmids) # calcuate total no of abstracts
        GeneObj.get_pmid_pos(annoDB=annDB)
        abs_in_DB += len(GeneObj.cellinpmid)
        occuDF = GeneObj.get_occurrence(cellDB=occuDF)
        count -= 1
    joincellsynonym(occuDF, resource_path)
    if not created:
        occuDF.to_csv(resource_path + os.path.sep + 'gene_occu_db.csv', sep=',', ignore_index=True)
    if join:
        print("\nUpdating gene occurrence db....")
        update_dataframe = pd.concat([occuDF.drop(['celltype'], axis=1), dataframe], axis=1)
        update_dataframe.to_csv(resource_path + os.path.sep + 'gene_occu_db.csv', sep=',', ignore_index=True)
    occuDF = pd.concat([occuDF, foundgenes_df], axis=1)
    return occuDF

def joincellsynonym(celloccu, resource_path):
    '''
    Join multiple cell synonym to one.
    :param celloccu:
    :param cellSyn:
    :return:
    '''
    cellSyn = FilesFolders.cell_synonym(resource_path)
    colname = celloccu.columns.values.tolist()
    indexRem = []
    #print celloccu
    for k, v in cellSyn.iterrows():
        index = celloccu.celltype[celloccu.celltype == v['cell'].lower()].index.tolist()[0]
        for cell in v['synonyms'].split(','):
            #print cell
            indexsyn = celloccu.celltype[celloccu.celltype == cell.lower()].index.tolist()[0]
            indexRem.append(indexsyn)
            ## adding synonym
            for col in colname:
                if col != 'celltype' and col != 'Unnamed: 0':
                    celloccu.loc[index, col] = celloccu.loc[index, col] + celloccu.loc[indexsyn, col]
    celloccu = celloccu.drop(celloccu.index[indexRem])
    #print celloccu
    return celloccu