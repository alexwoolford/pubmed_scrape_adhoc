
file_handle_articles = open('doctor_articles.txt', 'r')

for line in file_handle_articles.readlines():
    doctor_name, pmid, journal_title, article_title, abstract = line.strip().split('|')
    print doctor_name, pmid
