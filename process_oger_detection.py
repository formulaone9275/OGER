from __future__ import unicode_literals, print_function
import sys,codecs
import os,csv,re
from nltk.tokenize import word_tokenize
from nltk.corpus import words




def detect_terms_after_c_term_(text):
    list_of_word=text.split(' ')
    dict= words.words()

    word_extra=[]

    for wii in word_extra:
        dict.append(wii)

    #if there is a capital letter or a digit in the word, we will see it as c-term
    pattern_for_capital_letter='[A-Z]+'
    pattern_for_digit='\d+'
    pattern_for_any_letter='[A-Za-z]+'
    com_pattern_capital=re.compile(pattern_for_capital_letter)
    com_pattern_digit=re.compile(pattern_for_digit)
    com_pattern_any_letter=re.compile(pattern_for_any_letter)
    c_term_index=-1
    for i in range(len(list_of_word)):
        wi=list_of_word[i]
        wi=wi.replace('+','')
        wi=wi.replace('*','')
        if len(wi)==0:
            continue


        if wi.lower() not in dict:
            sr_capital=com_pattern_capital.search(wi)
            sr_digit=com_pattern_digit.search(wi)
            sr_letter=com_pattern_any_letter.search(wi)
            if sr_capital or (sr_digit and sr_letter):
                c_term_index=i
                break

    if c_term_index<=0:
        return '',text

    return ' '.join(list_of_word[:c_term_index]),' '.join(list_of_word[c_term_index:])


def run():
    csvfile='uniprot_human_sprot.csv'
    change_list=[]
    csvfile_new='c_term_protein_name.csv'

    word_before_c_term_pattern='(activator|\w*ase|\w*in|protein|.*antigen|subunit|adaptor|interactor|inhibitor' + \
                               '|activator|.*factor|component|chaperone|cytochrome|suppressor|.*receptor|histone|member)\s?$'

    c_term_start_pattern='^(\d+|[A-Za-z]|[XIV]+)(\W+|$)'
    complie_word_before_c_term=re.compile(word_before_c_term_pattern)
    complie_c_term_start=re.compile(c_term_start_pattern)

    with codecs.open(csvfile_new, 'w',encoding='utf-8') as tsvfile1:
        spamwriter = csv.writer(tsvfile1, delimiter='\t', quotechar='|')
        col_title=['original name','new name','matched f term','droped name','rules applied']
        spamwriter.writerow(col_title)
        with codecs.open(csvfile, 'r') as tsvfile2:
            spamreader = csv.reader(tsvfile2, delimiter='\t', quotechar='|')

            for row in spamreader:
                #[original name, new name, matched f_term before cterm, dropped ones, rule working(Yes/No)]
                res_row=['']*5
                res_row[-1]='No'
                original_protein_name=row[3]
                words_before_c_term,new_protein_name=detect_terms_after_c_term_(original_protein_name)

                if new_protein_name!=original_protein_name:

                    res_row[0]=original_protein_name
                    res_row[1]=new_protein_name
                    sr_f_term=complie_word_before_c_term.search(words_before_c_term)
                    sr_new_protein=complie_c_term_start.search(new_protein_name)
                    #first check we find some special cases of c term like straing with digit, single letter or roman letter
                    #we just write the results in the dropped cases and continue to next one
                    if sr_new_protein:
                        res_row[3]=new_protein_name
                        if res_row not in change_list:
                            change_list.append(res_row)
                            spamwriter.writerow(res_row)
                        continue

                    #then match the word before the c-term
                    if sr_f_term:
                        res_row[2]=words_before_c_term[sr_f_term.span()[0]:sr_f_term.span()[1]]
                        res_row[-1]='Yes'

                    else:
                        res_row[3]=new_protein_name
                        #res_row[-1]='Yes'

                    if res_row not in change_list:
                        change_list.append(res_row)
                        spamwriter.writerow(res_row)
                    #print(row)
                    print(res_row)



if __name__ == '__main__':
    text='Plasma protease C1 inhibitor'
    print(detect_terms_after_c_term_(text))
    run()