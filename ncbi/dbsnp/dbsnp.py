#!/usr/bin/env python3

import time
import ncbiutils
from Bio import Entrez


"""
Finds all rsids that are explicitly cited in pubmed
and returns a list
"""
def get_complete_rsids():
    results = ncbiutils.db_query(db="snp",retmode="json",retmax=200000,retstart=0,term='snp_pubmed_cited[sb]')
    rsidlist = results["esearchresult"]["idlist"]
    for x in range(0, len(rsidlist)):
        rsidlist[x] = "rs" + rsidlist[x]
    return(rsidlist)

"""
Generates a list of PMIDs that are explicitly cite a given rsid
"""
def get_pmids(rsid):
    searchterm = rsid + "+AND+pubmed_snp_cited[sb]"
    print(searchterm)
    results = ncbiutils.db_query(db="pubmed",retmode="json",retmax=200000,restart=0,term=searchterm, api_key="7c0213f7c513fa71fe2cb65b4dfefa76fb09")
    pmidlist = results["esearchresult"]["idlist"]
    print(pmidlist)
    return(pmidlist)


"""
Takes the saved list of results from get_pmids and retrieves
the raw XML to parse.  Currently goes article by article instead of
doing them in bulk.  Could go either way. Unsure which way
is better / more efficient.
"""

def get_abstracts(results):
    abstracts_list = []
    for start in range(0, int(results["Count"]), 1):
        # print("Going to download record %i to %i" % (start+1, end))
        fetch_handle = Entrez.efetch(db="pubmed", rettype="abstract",
                                        retmode="xml", retstart=start,
                                        retmax=1,
                                        webenv=results["WebEnv"],
                                        query_key=results["QueryKey"])
        data = fetch_handle.read()
        fetch_handle.close()
        root = ET.fromstring(data)
        for abst in root.iter('Abstract'):
            for sec in abst.iter('AbstractText'):
                abstracts_list.append(sec.text)
    print(abstracts_list)
    return(abstracts_list)


def get_abstracts_from_list(pmids_list):
    abstracts_list = []
    for each_pmid in pmids_list:
        fetch_handle = Entrez.efetch(db="pubmed", id=each_pmid, retmode='xml')
        data = fetch_handle.read()
        fetch_handle.close()
        root = ET.fromstring(data)
        for abst in root.iter('Abstract'):
            for sec in abst.iter('AbstractText'):
                abstracts_list.append(sec.text)
    return(abstracts_list)

def main():
    rsids = get_complete_rsids()
    for rsid in rsids:
        get_pmids(rsid)
    return()

if __name__ == '__main__':
    main()