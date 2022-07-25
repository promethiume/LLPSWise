# -*- coding: utf-8 -*-
# ************************************************************ #
#   FileName      : uniprot_function.py
#   Author        : Mengchen
#   Email         : mengchenpu@gmail.com
#   Create on     : 06-21-2022
#   Last modified : 06-21-2022 17:22
#   Version       : V1.0
#   Description   : 
# ************************************************************ #
import re
import time
import urllib
from Bio import SeqIO
from bs4 import BeautifulSoup

def get_function(query=""):
    """
    input: uniprotID
    output: function list
    
       get protein function from uniport (original from GO annotation)
    """
    function_list=[]
    go_celllocation_list = []
    go_pathway_list = []
    
    with open("uniprotDB2/"+query+".txt", "r") as f:
        page = f.readlines()

    for line in page:
        if 'DR   GO;' in line and '; F:' in line:
            function_list.append(line.split("; ")[2][2:])
        if 'DR   GO;' in line and '; C:' in line:
            go_celllocation_list.append(line.split("; ")[2][2:])
        if 'DR   GO;' in line and '; P:' in line: 
            go_pathway_list.append(line.split("; ")[2][2:])
 
    #print(function_list, go_celllocation_list, go_pathway_list)
    return function_list, go_celllocation_list, go_pathway_list


def uniprot(uniprotID):
    """
    input: uniprotID
    output: sequence, subcellular location (list), genename
 
        get subcellular location
    """
    for record in SeqIO.parse("uniprotDB2/"+uniprotID+".txt", "swiss"):
        resultseq = record.seq

        try:
            resultGeneName = re.split(" |{|;", record.annotations['gene_name'])[0].split("=")[1]
        except:
            resultGeneName = "-"

        try:
            rec = re.split(r"\n", record.annotations['comment'].split("SUBCELLULAR LOCATION:")[1].split("Note")[0])[0]
            resultLocation = re.sub('\{[a-zA-Z0-9:, |\-]+\}','',rec).lstrip().split(". ")
            for index, i in enumerate(resultLocation):
                i_new=re.sub('\.', '', i)
                resultLocation[index]=i_new.rstrip()
        except:
            resultLocation = []

        if len(resultLocation)>1 and len(resultLocation[-1])==0:
            resultLocation.pop()

    return resultseq, resultLocation, resultGeneName

#a,b,c=uniprot("Q9UER7")
#print(a,b,c)
