__author__ = 'peeyush'


class Genes():
    def __init__(self, gene, subquery=None):
        self.gene = gene
        self.subquery = subquery
        self.pmids = None
        self.cellinpmid = None
        self.occurrence = None


    def get_pmids(self):
        from Bio import Entrez
        gene = self.gene.strip()
        subquery = self.subquery
        Entrez.email = "peeyush215@gmail.com"
        if subquery:
            newsubquery = subquery.strip().replace(',', ' AND ')
            query = gene +' AND '+ newsubquery
            data = Entrez.esearch(db="pubmed", retmax=15000, term=query)
            res = Entrez.read(data)
            PMID = res["IdList"]
            #print 'Pubmed ids for '+query+':', len(PMID)
            self.pmids = PMID
        else:
            data = Entrez.esearch(db="pubmed", retmax=15000, term=gene)
            res = Entrez.read(data)
            PMID = res["IdList"]
            #print 'Pubmed ids for '+gene+':', len(PMID)
            self.pmids = PMID


    def get_pmid_pos(self, annoDB):
        celltype_list = []
        for pmid in self.pmids:
            if int(pmid) in annoDB.index:
                celltype_list.append(annoDB.loc[int(pmid)][0])
        self.cellinpmid = celltype_list


    def get_occurrence(self, cellDB):
        '''
        Calculate celltype occurrence for each gene.
        :param genes_dict:
        :param cellDB:
        :return:
        '''
        celloccu = cellDB
        # print cellDB
        celloccu[self.gene] = 0
        #print 'Gene:', self.gene
        celltype = self.cellinpmid
        for found in celltype:
            for index, cells in celloccu.iterrows():
                if cells['celltype'].strip().lower() in found.lower():
                    celloccu.loc[index, self.gene] += 1
        celloccu['celltype'] = celloccu['celltype'].str.lower()
        return celloccu

def pmid_binary_search(annoDB, pmid, imin, imax):
    '''
    This function will search the index of pmid.
    Caution: To perform binary search you always need a index sorted dataframe.
    :param annoDB:
    :param pmid:
    :param imin:
    :param imax:
    :return:
    '''
    if imax < imin:
        return 'PMID NOT FOUND:'+str(pmid)
    else:
        imid = imin + ((imax - imin)/2)
        if annoDB.iloc[imid][0] > pmid:
            return pmid_binary_search(annoDB, pmid, imin, imid-1)
        elif annoDB.iloc[imid][0] < pmid:
            return pmid_binary_search(annoDB, pmid, imid+1, imax)
        else:
            return imid
