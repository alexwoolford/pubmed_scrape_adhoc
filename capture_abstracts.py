#!/usr/bin/env python

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


def capture_abstracts():
    doctor_names_filehandle = codecs.open("doctor_names.txt", "r", "utf-8")
    doctor_articles_filehandle = codecs.open("doctor_articles.txt", "w", "utf-8")

    for doctor_name in doctor_names_filehandle.readlines():

        doctor_name = doctor_name.strip()
        time.sleep(1)
        id_list = search(doctor_name)

        if id_list:

            id_details = fetch_details(id_list)
            if id_details:

                pubmed_articles = id_details['PubmedArticle']
                for pubmed_article in pubmed_articles:

                    try:
                        pmid = pubmed_article['MedlineCitation']['PMID']
                    except:
                        pmid = ""
                    try:
                        article_title = pubmed_article['MedlineCitation']['Article']['ArticleTitle']
                    except:
                        article_title = ""
                    try:
                        journal_title = pubmed_article['MedlineCitation']['Article']['Journal']['Title']
                    except:
                        journal_title = ""
                    try:
                        article_date_dict = pubmed_article['MedlineCitation']['Article']['ArticleDate'][0]
                        article_date = "{0}/{1}/{2}".format(article_date_dict.get('Month'), article_date_dict.get('Day'), article_date_dict.get('Year'))
                    except:
                        article_date = ""
                    try:
                        abstract = " ".join(pubmed_article['MedlineCitation']['Article']['Abstract']['AbstractText'])
                    except:
                        abstract = ""

                    doctor_articles_filehandle.write("|".join([doctor_name, pmid, journal_title, article_title, article_date, abstract]) + "\n")


if __name__ == "__main__":
    capture_abstracts()
