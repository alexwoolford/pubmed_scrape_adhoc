#!/usr/bin/env python

import excel_doctor_extract
from Bio import Entrez
import time
import codecs


def search(query):
    try:
        Entrez.email = 'alex@woolford.io'
        handle = Entrez.esearch(db='pubmed',
                                sort='relevance',
                                retmax='200',
                                retmode='xml',
                                term=query)
        results = Entrez.read(handle)
        if results['IdList']:
            return results['IdList']
        else:
            return None
    except:
        return None


def fetch_details(id_list):
    try:
        ids = ",".join(id_list)
        Entrez.email = 'alex@woolford.io'
        handle = Entrez.efetch(db='pubmed',
                               retmode='xml',
                               id=ids)
        results = Entrez.read(handle)
        return results
    except:
        return None


if __name__ == "__main__":

    with codecs.open("doctor_articles.txt", "w", "utf-8") as file_handle:

        column_names = ['doctor_name', 'pmid', 'journal_title', 'article_title', 'abstract']

        file_handle.write("|".join(column_names) + "\n")

        excel_doctor_extract = excel_doctor_extract.ExcelDoctorExtract()
        doctor_names = excel_doctor_extract.parse('Faculty TOBI 2010-2016 & 2017.xlsx')
        for doctor_name in doctor_names:
            time.sleep(1)
            id_list = search(doctor_name)
            if id_list:
                id_details = fetch_details(id_list)
                if id_details:
                    for id_detail in id_details:
                        try:
                            pmid = id_detail['MedlineCitation']['PMID']
                        except:
                            pmid = ""
                        try:
                            article_title = id_detail['MedlineCitation']['Article']['ArticleTitle']
                        except:
                            article_title = ""
                        try:
                            journal_title = id_detail['MedlineCitation']['Article']['Journal']['Title']
                        except:
                            journal_title = ""
                        try:
                            abstract = " ".join(id_detail['MedlineCitation']['Article']['Abstract']['AbstractText'])
                        except:
                            abstract = ""
                        file_handle.write("|".join([doctor_name, pmid, journal_title, article_title, abstract]) + "\n")

