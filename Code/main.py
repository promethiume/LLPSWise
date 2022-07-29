#   FileName      : main-origllps.py
#   Author        : Mengchen
#   Email         : mengchenpu@gmail.com
#   Create on     : 02-25-2022
#   Last modified : 03-17-2022 14:18
#   Version       : V1.0
#   Description   : inherited from data_extractionV4.py 
# ************************************************************ #

import warnings
warnings.filterwarnings("ignore")

import esm
import ssl
import sys
import pickle
import urllib
import argparse
import subprocess

import numpy as np
import pandas as pd
import networkx as nx

from Bio import SeqIO
from bs4 import BeautifulSoup

import torch
import torchvision
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import Variable
from torch.utils.data import Dataset

from itertools import chain  
import matplotlib.pyplot as plt
from pyecharts import options as opts            
from pyecharts.charts import Graph as Graph

from model import ANNmodel
from uniprot_function3 import get_function, uniprot 

parser = argparse.ArgumentParser(description='Process input protein ID (uniprotID)')
parser.add_argument('--id', help='input uniprotID')
uniprotID = parser.parse_args().id

# set a couple of parameters
priConN=10
secConN=5
scaffold_priConN = 10
scaffold_secConN = 5
llpsCutoff = 0.85
bootstrap = "union"   # union/no
bootstrapN = 3
maxEdges = 4
binding=[]
blackList = ["P0CG48"] #UBC, 

print(bootstrap)

# set up protein-protein interaction database
ppidata = pd.read_table("BIOGRID-ORGANISM-Homo_sapiens-4.4.207.tab3.txt", low_memory=False)

flatten = lambda x: [i for row in x for i in row] # flatten 2D list to 1D

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
 
input_dim = 1280
hidden_dim = 64                                                                 
output_dim = 2
model = ANNmodel(input_dim, hidden_dim, output_dim)
model=torch.load("model.pth")
model.to(device)
model.eval()

# Load ESM-1b model and get sequence vector
seqmodel, alphabet = esm.pretrained.esm1_t34_670M_UR50S()
batch_converter = alphabet.get_batch_converter()
 
# init nn for phase separation prediction from target seq
def llpsp(targetseq):
    if len(targetseq) <= 2000:
        seqData =[("targetSequence", targetseq)]
        batch_labels, batch_strs, batch_tokens = batch_converter(seqData)
        with torch.no_grad():
            results = seqmodel(batch_tokens, repr_layers=[34], return_contacts=True)
        token_representations = results["representations"][34]
     
        for j, (_, seq) in enumerate(seqData):
            seqVector = token_representations[j, 1 : len(seq) + 1].mean(0)
        
        inputfile = seqVector.unsqueeze(0).to(device)
        output = model(inputfile)
        accuracy = (output.argmax(1)==0)
        prob = F.softmax(output, dim=1).cpu().tolist()[0][0]
    else:
        prob = 0.0

    print(prob)

    return prob, "pass"


def get_ppidata(targetGene, ID, parN, ResultTable, pri="yes"):
    """
    input: geneName, uniprotID, ppi-circle-size, calculate llps(Y/N)
    output: contact gene name(list), contact uniprotID(list), 
            if calculate llps(calculate primaray contact)
    """
    # get 1st contact protein
    _1contactID=ppidata[(ppidata["SWISS-PROT Accessions Interactor A"]==ID)& (ppidata['Organism Name Interactor B']=="Homo sapiens")][["SWISS-PROT Accessions Interactor B", "Publication Source", 'Throughput', 'Experimental System', 'Official Symbol Interactor B']]
    _2contactID=ppidata[(ppidata["SWISS-PROT Accessions Interactor B"]==ID)& (ppidata['Organism Name Interactor A']=="Homo sapiens")][["SWISS-PROT Accessions Interactor A", "Publication Source", 'Throughput', 'Experimental System', 'Official Symbol Interactor A']]
    _1contactID.columns = ["SWISS-PROT Accessions Interactor A", "Publication Source", 'Throughput', 'Experimental System', 'Official Symbol Interactor A']
    _2contactID.columns = ["SWISS-PROT Accessions Interactor A", "Publication Source", 'Throughput', 'Experimental System', 'Official Symbol Interactor A']
    contactID=pd.concat([_1contactID, _2contactID], axis=0)
    contactID.drop_duplicates(keep="first", inplace=True)
    print(contactID[["SWISS-PROT Accessions Interactor A", 'Official Symbol Interactor A']].value_counts().head(parN))
    contactID=contactID[["SWISS-PROT Accessions Interactor A", 'Official Symbol Interactor A']].value_counts().head(parN).index
    contactID=contactID.to_frame(index=False)['SWISS-PROT Accessions Interactor A'].tolist()

    try:
        contactID.remove("-")
    except:
        pass
    try:
    #    contact.remove(targetGene)
        contactID.remove(ID)
    except:
        pass
    
    #print(contact)
    resContact=[]
    resContactID=[]
    pcllpsScore=0.0
    pcllpsTop=""
    pcllpsTopID=""
    ContactScaffold=[]
    ContactScaffoldID=[]
    llpsList=[]
    keepFlag = "-"
    bindingFlag = "-"
    

    for i, iID in enumerate(contactID):
        iID = iID.split("|")[0]

        if ResultTable[ResultTable['uniprotID']==iID].empty:
            pcseq, pcLocation, pcGeneName = uniprot(iID)
            pcFunction, pcGoLocation, pcGoPathway = get_function(iID)
            llpspCanScore, _ = llpsp(pcseq)
            print('NEEEEEEW!!!!')
        else:
            pcseq = ResultTable[ResultTable['uniprotID']==iID]['seq'].item()
            pcLocation = ResultTable[ResultTable['uniprotID']==iID]['location'].item()
            pcGeneName = ResultTable[ResultTable['uniprotID']==iID]['geneName'].item()
            pcFunction = ResultTable[ResultTable['uniprotID']==iID]['function'].item()
            pcGoLocation = ResultTable[ResultTable['uniprotID']==iID]['golocation'].item()
            pcGoPathway = ResultTable[ResultTable['uniprotID']==iID]['goPathway'].item()
            llpspCanScore = ResultTable[ResultTable['uniprotID']==iID]['llpsScore'].item()
        print(pcGeneName, iID, pcLocation, pcFunction)

        if (list(set(pcLocation)&set(targetLocation)) or list(set(pcGoLocation)&set(targetGoLocation))) and (list(set(pcFunction)&set(targetFunction)) or list(set(pcGoPathway)&set(targetGoPathway)) ):
            print("YES, share location and function")
            resContact.append(pcGeneName)
            resContactID.append(iID)

            keepFlag = "Y"
            
            llpsList.append([iID, llpspCanScore])
            if "RNA binding" in str(pcFunction) or "DNA binding" in str(pcFunction) or "RNA binding" in str(pcGoPathway) or "DNA binding" in str(pcGoPathway):
                binding.append(pcGeneName)
                bindingFlag = "Y"
            else:
                bindingFlag = "N"
            
            if pri=="yes":
                print(llpspCanScore)
                if float(llpspCanScore) > float(pcllpsScore):
                    pcllpsScore = float(llpspCanScore)
                    pcllpsTop, pcllpsTopID = pcGeneName, iID
                else:
                    pass
                ContactScaffold.append(pcGeneName)
                ContactScaffoldID.append(iID)
                #if float(llpspCanScore) > float(llpsCutoff):
                #    ContactScaffold.append(pcGeneName)
                #    ContactScaffoldID.append(contactID[i])
            else:
                pcllpsScore, pcllpsTop, pcllpsTopID ="NA", "NA","NA"
            
        else:
            print(pcLocation, targetLocation)
            print(pcGoLocation, targetGoLocation)
            print(pcFunction, targetFunction)
            print(pcGoPathway, targetGoPathway)
            print("NOT share function and location!")
            keepFlag = "N"
        
        if ResultTable[ResultTable['uniprotID']==iID].empty:
            size = ResultTable.index.size
            #print(size)
            ResultTable.loc[size] = [pcGeneName, iID, pcseq, llpspCanScore, pcLocation, pcGoLocation, pcFunction, pcGoPathway, keepFlag, bindingFlag]
            #print(ResultTable)

                  
    print(resContact, resContactID)
    return resContact, resContactID, pcllpsTop, pcllpsTopID, pcllpsScore, ContactScaffold, ContactScaffoldID, llpsList, ResultTable

queryResult = pd.DataFrame(columns = ['geneName', 'uniprotID', 'seq', 'llpsScore', 'location', 'golocation', 'function', 'goPathway', 'shareTS', 'nucbinding'])

iseq, targetLocation, geneName = uniprot(uniprotID)
print(geneName)
targetFunction, targetGoLocation, targetGoPathway = get_function(uniprotID)
llpspScore, llpsYN = llpsp(iseq)
queryResult.loc[0] = [geneName, uniprotID, iseq, llpspScore, targetLocation, targetGoLocation, targetFunction, targetGoPathway, 'target', 'target']
#print(queryResult)

priContact, priContactID, llpsTop, llpsTopID, topScore, priScaffold, priScaffoldID, prillpsList, queryResult = get_ppidata(geneName, uniprotID, priConN, queryResult)
#print(priContact)
#print(queryResult)

mainGraph = nx.Graph()

print(llpspScore)

if float(llpspScore) > float(llpsCutoff):  # target protein is scaffold
    print(geneName, " propensity: ", llpsYN, llpspScore)
    print(targetLocation, targetFunction)
    print("priContact:", priContact)
    nodelist1 = [geneName]
    graphInputEdge = [[a, b] for a in nodelist1 for b in priContact if a != b]
    mainGraph.add_edges_from(graphInputEdge)

    

    secContact=[]
    secllpsList=[]
    for i, iID in enumerate(priScaffold):
        iContact, _, _, _, _ , _ ,_,isecllpsList, queryResult= get_ppidata(priContact[i], priContactID[i], secConN, queryResult, "No")
        secllpsList = secllpsList+isecllpsList
        secContact.append(iContact)
        nodelist1 = [priContact[i]]              
        graphInputEdge = [[a, b] for a in nodelist1 for b in iContact if a != b]
        mainGraph.add_edges_from(graphInputEdge)
    print("secContact:", secContact)
    flag="Yes"

elif float(topScore) > float(llpsCutoff):   # target protein is client but pri contact include a scaffold protein
    print(geneName, " propensity: ", llpsYN, llpspScore)
    print(targetLocation, targetFunction)
    nodelist1 = [geneName]              
    graphInputEdge = [[a, b] for a in nodelist1 for b in priContact if a != b]
    mainGraph.add_edges_from(graphInputEdge)

    priContact, priContactID, _,_,_, priScaffold, priScaffoldID, prillpsList, queryResult = get_ppidata(llpsTop, llpsTopID, priConN, queryResult)
    print(llpsTop, " propensity: ", topScore)
    print("priContact:", priContact)
    nodelist1 = [llpsTop]              
    graphInputEdge = [[a, b] for a in nodelist1 for b in priContact if a != b]
    mainGraph.add_edges_from(graphInputEdge)

    secContact=[]     
    secllpsList=[]
    for i, iID in enumerate(priScaffold): 
        iContact, _, _,_,_, _, _, isecllpsList, queryResult = get_ppidata(priContact[i], priContactID[i], secConN, queryResult, "No")
        secllpsList = secllpsList+isecllpsList
        secContact.append(iContact)
        nodelist1 = [priContact[i]]                                                                              
        graphInputEdge = [[a, b] for a in nodelist1 for b in iContact if a != b]                                 
        mainGraph.add_edges_from(graphInputEdge)
    print("secContact:", secContact)
    flag="Help"
else:    # target protein is client and pri contact DONOT include a scaffold protein
    print("phase separation protensity cut off:", llpsCutoff)
    print(geneName, " propensity: ", llpspScore)
    print(targetLocation, targetFunction)
    print(llpsTop, " propensity: ", topScore)
    print("No")
    c = (     
    Graph()                   
    .set_global_opts(title_opts=opts.TitleOpts(title="low bio-condensate propensity"))  
    .render("templates/newnetworkpriC3-0616"+str(uniprotID)+bootstrap+".html")
    )
    print("Not a LLPS system")
    sys.exit()

# node generator
target = [geneName]
if flag == "Help":
    scaffold = [llpsTop]
else:
    scaffold = []

# start bootstrapping!
bsGraph = nx.Graph()

scaffoldList=[]
if bootstrap != "no":
    llps_propensity_list=[]
    llps_propensity_list.append([uniprotID, llpspScore])
    print(prillpsList)
    print(secllpsList)
    llps_propensity_list = llps_propensity_list + prillpsList + secllpsList
    print(llps_propensity_list)
    Topscaffold = pd.DataFrame(llps_propensity_list).drop_duplicates()
    
    Topscaffold = Topscaffold[~Topscaffold[0].isin(blackList)]

    Topscaffold.sort_values(by=1, inplace=True, ascending=False)
    print(Topscaffold)
    print(Topscaffold.head(bootstrapN))

    orig_nodes=list(set(target+scaffold+priContact+flatten(secContact)))
    
    scaffold_priContact, scaffold_priContactID, scaffold_secContact=[], [], []
    
    
    driverList = [geneName]
    for index, i in enumerate(Topscaffold.head(bootstrapN)[0]):
        if i != uniprotID and float(Topscaffold.head(bootstrapN)[1].iloc[index]) > llpsCutoff:
            print("!!!start bootstrapping")
            print(Topscaffold.head(bootstrapN)[1].iloc[index])
    
            iscaffoldSeq, iscaffoldLocation, iscaffold_geneName = uniprot(i)
            scaffoldList.append(iscaffold_geneName)
            iscaffold_Function, iscaffold_GoLocation, iscaffold_GoPathway = get_function(i)
           
            iscaffold_priContact, iscaffold_priContactID, _, _, _, _, _ , _, queryResult = get_ppidata(iscaffold_geneName, i, scaffold_priConN, queryResult, "Yes")

            nodelist1 = [iscaffold_geneName]                    
            graphInputEdge = [[a, b] for a in nodelist1 for b in iscaffold_priContact if a != b]
            bsGraph.add_edges_from(graphInputEdge)
            driverList.append(iscaffold_geneName)

            if flag == "Help":
                driverList.append(llpsTop)
            
            iscaffoldSecContact=[]
            iscaffold_priContact_inter=[]
            iscaffold_priContactID_inter=[]
            for j, jID in enumerate(iscaffold_priContactID):
                jContact, _, _, _, _ , _ ,_, _, queryResult= get_ppidata(iscaffold_priContact[j], iscaffold_priContactID[j], scaffold_secConN, queryResult, "No")
                iscaffoldSecContact.append(jContact)

                nodelist1 = [iscaffold_geneName]                                  
                graphInputEdge = [[a, b] for a in nodelist1 for b in iscaffold_priContact if a != b]
                bsGraph.add_edges_from(graphInputEdge)

    unionGraph = nx.compose(bsGraph, mainGraph)
    print(unionGraph.nodes)
    print(driverList)
    priConList = []
    for inodes in driverList:   # generate priContact list
        path = nx.single_source_shortest_path(unionGraph, inodes)
        for ikey in path:
            if len(path[ikey])==2:
                priConList.append(ikey) 
    
    for inodes in driverList:
        try:
            path = nx.single_source_shortest_path(unionGraph, inodes)
            print(inodes)
            for ikey in path:
                try:
                    if len(path[ikey])> maxEdges and ikey not in priConList and ikey not in driverList:
                        print(path[ikey])
                        print(ikey)
                        unionGraph.remove_node(ikey)
                except:
                    pass
        except:
            pass

print("generate graph")

if bootstrap != "no":
    contact = unionGraph.nodes
else:
    contact = mainGraph.nodes

binding = list(set(binding))
nodes, links=[], []  
for i in target:              
    nodes.append({"name": i, "symbolSize": 40, "itemStyle": {"normal": {"color": "#FFA500", "borderColor": "red", "borderWidth": 5}}})
          
if scaffold:                  
    for j in scaffold:  
        nodes.append({"name": j, "symbolSize": 35, "itemStyle": {"normal": {"color": "orange"}}})

if scaffoldList:
    for j in scaffoldList:
        if j not in target and j not in scaffold: 
            nodes.append({"name": j, "symbolSize": 30, "itemStyle": {"normal": {"color": "#FF6347"}}})

if binding:
    for i in binding:
        if i not in scaffold and i not in scaffoldList and i not in target and i in contact:
            nodes.append({"name": i, "symbolSize": 20, "itemStyle": {"normal": {"color": "#1E90FF", "borderColor": "#1E9000", "borderWidth": 5}}})

if bootstrap == "no":
    for i in mainGraph.nodes:       
        if i not in target and i not in scaffold and i not in scaffoldList and i not in binding:                   
            nodes.append({"name": i, "symbolSize": 20, "itemStyle": {"normal": {"color": "#1E90FF"}}})

if bootstrap != "no":
    for i in unionGraph.nodes:       
        if i not in target and i not in scaffold and i not in scaffoldList and i not in binding:                   
            nodes.append({"name": i, "symbolSize": 20, "itemStyle": {"normal": {"color": "#1E90FF"}}})


# link generator
if bootstrap == "no":
    for i in mainGraph.edges:
        print(i)
        links.append({"source": i[0], "target": i[1]})

if bootstrap != "no":
    for i in unionGraph.edges:
        print(i)
        links.append({"source": i[0], "target": i[1]})

c = (
    Graph()
    .add("", nodes, links, repulsion=4000)
    .set_global_opts(title_opts=opts.TitleOpts(title="Bio-condensate Network of "+geneName+" "+bootstrap ))  
    #.render("templates/graph_with_options.html") 
    .render("templates/newnetworkpriC3-0616"+str(uniprotID)+bootstrap+".html") 

)
queryResult.to_csv(str(uniprotID)+"-"+str(geneName)+".csv")

