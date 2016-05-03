__author__ = 'peeyush'
from GCAM import FilesFolders

def gene2synonym(geneList, geneSyn):
    '''
    This will retrieve synonyms for user provided genes
    :param geneList:
    :param geneSyn:
    :return:
    '''
    newGeneList = []
    for gene in geneList:
        newGeneList.append(gene)
        if gene in geneSyn.index:
            synonym = geneSyn.loc[gene][1].split(',')
            for syn in synonym:
                if len(syn.strip()) > 0 and str(syn.strip()) != '0':
                    newGeneList.append(syn.strip().lower())
    return newGeneList


def get_occurrence(genes_dict, cellDB):
    '''
    Calculate celltype occurrence for each gene.
    :param genes_dict:
    :param cellDB:
    :return:
    '''
    celloccu = cellDB
    for k, v in genes_dict.iteritems():
        celloccu[k] = 0
        celltype = v.cellinpmid
        for found in celltype:
            for index, cells in celloccu.iterrows():
                if cells['celltype'] in found.lower(): #if cells['celltype'].strip().lower() in found.lower():
                    celloccu.loc[index, k] += 1
    celloccu['celltype'] = celloccu['celltype'].str.lower()
    return celloccu

def joincellsynonym(celloccu, cellSyn):
    '''
    Join multiple cell synonym to one.
    :param celloccu:
    :param cellSyn:
    :return:
    '''
    colname = celloccu.columns.values.tolist()
    indexRem = []
    #print colname
    for k, v in cellSyn.iterrows():
        index = celloccu.celltype[celloccu.celltype == v['cell'].lower()].index.tolist()[0]
        for cell in v['synonyms'].split(','):
            #print cell
            indexsyn = celloccu.celltype[celloccu.celltype == cell.lower()].index.tolist()[0]
            indexRem.append(indexsyn)
            for col in colname:
                if col != 'celltype' and col != 'Unnamed: 0':
                    celloccu.loc[index, col] = celloccu.loc[index, col] + celloccu.loc[indexsyn, col]
    celloccu = celloccu.drop(celloccu.index[indexRem])
    return celloccu

def joingenesynonym(colloccu, primarygemename, geneSyn):
    '''
    Join gene synonyms to one
    :param colloccu:
    :return:
    '''
    #print 'Shape of df before gene merge:', colloccu.shape
    col2drop = []
    for gene in primarygemename:
        if gene in geneSyn.index and geneSyn.loc[gene][1] != '0':
            synonym = geneSyn.loc[gene][1].strip(',').split(',')
            for syn in synonym:
                #print colloccu[syn.lower()].shape
                if len(colloccu[syn.lower()].shape) > 1:
                    raise ValueError("Synonym Error: Gene list contains more than 2 gene synonym belogns to same gene. You can re-run without synonym parameter.")
                col2drop.append(syn.lower())
                colloccu[gene] = colloccu[gene] + colloccu[syn.lower()]
    colloccu = colloccu.drop(col2drop, axis=1)
    #print 'Shape of df after gene merge:', colloccu.shape
    return colloccu

def subtract_cellnamepeat(celloccu, path):
    '''
    This will subtract count from ex. t lymphocyte that belongs to cd4 t cells
    :return:
    '''
    subtract_df = FilesFolders.read_cell_subtractdf(path)
    colnames = celloccu.columns
    #print celloccu.index
    #print subtract_df
    for k, v in subtract_df.iterrows():
        for cell in v['subs'].split(','):
            cell = cell.lower()
            for gene in colnames:
                #print gene
                if celloccu.loc[v['celltype'], gene] > 0:
                    celloccu.loc[v['celltype'], gene] = celloccu.loc[v['celltype'], gene] - celloccu.loc[cell, gene]
    #print celloccu
    return celloccu