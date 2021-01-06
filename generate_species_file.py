from __future__ import division,print_function
import time,re,sys,os,ast,json
import pandas as pd
import pymongo
from collections import OrderedDict
from pymongo import MongoClient
from bson.son import SON
from bson.codec_options import CodecOptions
import re,csv,enchant


class monogo_db:
    def __init__(self,db,ent,text):

        self.db=db

        self.ent=ent

        self.text=text
        #for checking c-term

        #--- create database instances---
        # Environment variables
        mongodb_host = os.environ.get("MONGODB_HOST","0.0.0.0") # change to biotm2.cis.udel.edu before dockerizing
        mongodb_port = os.environ.get("MONGODB_PORT","27017")
        db_name = os.environ.get("DBNAME",self.db) # change database name for your own dbName
        textCollectionName = os.environ.get("COLLECTION_TEXT",self.text)
        entityCollectionName = os.environ.get("COLLECTION_ENTITY",self.ent)

        # Database URI
        MONGODB_URI = 'mongodb://'+mongodb_host+':'+mongodb_port+'/'

        # Database object
        client = MongoClient(MONGODB_URI)
        opts = CodecOptions(document_class=SON)

        # Database
        dbName = client[db_name] # glyco/unicarb

        # Collection
        self.textCollection = dbName[textCollectionName].with_options(codec_options=opts)
        self.entityCollection = dbName[entityCollectionName].with_options(codec_options=opts)


    def get_entity_list_from_mongoDB(self,pmid):


        protein_list=[]
        species_list=[]

        entityDoc = self.entityCollection.find_one({"docId":str(pmid)})
        if entityDoc:
            #print(entityDoc["entity"])


            for duid,entity in entityDoc["entity"].items():
                if entity["entityType"] == "Gene":
                    if len(entity["entityId"])>0:
                        protein_list.append((entity["entityText"],entity["charStart"],entity["charEnd"],entity["source"],entity["entityId"][0]["idString"]))
                    else:
                        protein_list.append((entity["entityText"],entity["charStart"],entity["charEnd"],entity["source"],''))
                if entity["entityType"] == "Protein":
                    if len(entity["entityId"])>0:
                        protein_list.append((entity["entityText"],entity["charStart"],entity["charEnd"],entity["source"],entity["entityId"][0]["idString"]))
                    else:
                        protein_list.append((entity["entityText"],entity["charStart"],entity["charEnd"],entity["source"],''))
                if entity["entityType"] == "Species":
                    if len(entity["entityId"])>0:
                        species_list.append((entity["entityText"],entity["charStart"],entity["charEnd"],entity["source"],entity["entityId"][0]["idString"]))
                    else:
                        species_list.append((entity["entityText"],entity["charStart"],entity["charEnd"],entity["source"],''))



        return protein_list,species_list




if __name__ == '__main__':
    #test

    db_class=monogo_db("uniprot","entities","text")

    pmidFile='uniprot.txt'
    pmidList = pd.read_csv(pmidFile,header=None).iloc[:,0].tolist()

    #pmidList=['9689040']
    pmidList=list(set(pmidList))
    count=0
    output_file='protein_with_species.csv'
    with open(output_file, 'w') as tsvfile2:
        spamwriter = csv.writer(tsvfile2, delimiter='\t', quotechar='|')
        spamwriter.writerow(['pmid','entity','entity_type','source','id'])
        for pi in pmidList:
            count+=1
            print(count,'-th:',pi)
            pl,sl=db_class.get_entity_list_from_mongoDB(pi)
            el=sl+pl
            print(el)
            entities=set()
            spamwriter.writerow([str(pi),'','','',''])
            for ei in sl:
                entities.add(('',ei[0],'Species',ei[3],ei[4]))
            for ei in pl:
                entities.add(('',ei[0],'Protein',ei[3],ei[4]))
            for line in entities:
                line=list(line)
                line=[li.encode('utf-8') for li in line ]
                spamwriter.writerow(list(line))