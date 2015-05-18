__author__ = 'peeyush'


def gene2synonym(geneList, geneSyn):
    '''
    This will retrieve synonyms for user provided genes
    :param geneList:
    :param geneSyn:
    :return:
    '''
    newGeneList = []
    geneSyn['gene'] = geneSyn['gene'].str.lower()
    geneSyn = geneSyn.set_index(geneSyn['gene'])
    for gene in geneList:
        #print gene
        newGeneList.append(gene)
        if gene in geneSyn.index:
            synonym = geneSyn.loc[gene][1].split(',')
            for syn in synonym:
                if len(syn.strip()) > 0 and str(syn.strip()) != '0':
                    newGeneList.append(syn.strip().lower())
    return newGeneList


def get_occurrence(genes_dict, cellDB):
    #Calculate celltype occurrence for each gene.
    celloccu = cellDB
    # print cellDB
    for k, v in genes_dict.iteritems():
        celloccu[k] = 0
        print 'Gene:', k
        celltype = v.cellinpmid
        for found in celltype:
            for index, cells in celloccu.iterrows():
                if cells['celltype'].strip().lower() in found.lower():
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
        index = celloccu.celltype[celloccu.celltype == v['cell'].strip().lower()].index.tolist()[0]
        for cell in v['synonyms'].split(','):
            indexsyn = celloccu.celltype[celloccu.celltype == cell.strip().lower()].index.tolist()[0]
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
    print 'Shape of df before gene merge:', colloccu.shape
    col2drop = []
    geneSyn['gene'] = geneSyn['gene'].str.lower()
    geneSyn = geneSyn.set_index(geneSyn['gene'])
    for gene in primarygemename:
        synonym = geneSyn.loc[gene][1].strip(',').split(',')
        for syn in synonym:
            col2drop.append(syn.lower())
            colloccu[gene] = colloccu[gene] + colloccu[syn.lower()]
    colloccu = colloccu.drop(col2drop, axis=1)
    print 'Shape of df after gene merge:', colloccu.shape
    return colloccu

