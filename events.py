# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 11:04:15 2016

@author: Steven Fletcher
"""

import re
import sys
import os

#function reads the sequences including ancestral ones
#INPUT - name + path of the seqfile with ancestors and there respective nodes 
#OUTPUT - list of tuples with (seqname,seq) 
def read_seqs(seqfile):
    f=open(seqfile,'r')
    lines=f.readlines()
    f.close()
    seqs={}
    for i in range(0,len(lines),2) : 
        if lines[i][0]=='>':
            seqs[lines[i].split('\n')[0]]=lines[i+1].split('\n')[0]   
    return seqs
    

    
#function reads the ancestral tree file format from the FastML results
#INPUT - name + path of the ancestral tree format file with 
#Ancestor file format - [name of the node, parent of the node, child/children of the node (if none empty column)]
#OUTPUT - list of the lines [name of the node, parent of the node, child/children of the node (if none empty column)] 
def read_tree_ancestorformat(ancestor_treefile):
    f=open(ancestor_treefile,'r')
    ancestortree=[]
    lines=f.read().splitlines()
    for line in lines :
        if not line or line[0]!='#':
            ancestortree.append(line.split())
    f.close()
    if not ancestortree[-1]:
        ancestortree.pop(0)
        ancestortree.pop()
    return ancestortree



#function finds the events(mutations) by comparing node with parent node
#INPUT - the sequences read from read_seqs(), the tree read for read_tree_ancestorformat()
#OUTPUT - a dictionary of events with the keys being the node name and the value 
#being the list of events given from the function pairwise_comparaison()
def find_events(seqs,ancestortree,i,AAlist):
   events={}
   for branch in ancestortree :
       if branch[1]!='root!' and i==0:
           events[branch[0]]=pairwise_comparaison(seqs['>'+branch[0]],seqs['>'+branch[1]])
       elif branch[1]!='root!' and i==1:
           events[branch[0]]=pairwise_custom_comparison(seqs['>'+branch[0]],seqs['>'+branch[1]],AAlist)
   return events



#function does the pairwise comparison between the sequence of the node and the parent node
#INPUT - the sequences of the node and the parent node
#OUTPUT - a list of events with 1 meaning there is a difference between the sequences and 0 none      
def pairwise_comparaison(seq_node,seq_parentnode):
    events=[]
    for pos in range(len(seq_node)):
        if seq_node[pos]!=seq_parentnode[pos] and seq_node[pos]!='-' and seq_parentnode[pos]!='-':
            events.append(1)
        else :
            events.append(0)
    return events



#function does the pairwise comparison between the sequence of the node and the parent node bit looks to see if the AAs are in a list
#in order to do custom events
#INPUT - the sequences of the node, the parent node and the AAlist
#OUTPUT - a list of events with 1 meaning there is a difference between the sequences and 0 none  
def pairwise_custom_comparison(seq_node,seq_parentnode,AAlist):
    events=[]
    for pos in range(len(seq_node)):
        if seq_node[pos] not in AAlist and seq_parentnode[pos] not in AAlist:
            events.append(0)
        elif seq_node[pos]==seq_parentnode[pos] :
            events.append(0)
        else :
            events.append(1)
    return events
 
   

#function creates the epics newick tree with the events 
#INPUT - the newick tree and the events
#OUTPUT - the epics newick tree      
def make_epicsNewick(Newick_tree_file,events):
    f=open(Newick_tree_file,'r')
    lines = f.read().splitlines() 
    f.close()
    epics_newick=lines[0]
    for key in events:
        if key+':' not in epics_newick:
            print 'skipped'
            continue
        else :
            pos=epics_newick.index(key+':')+len(key)
            event=format_event(events[key])
            epics_newick=epics_newick[:pos]+event+epics_newick[pos:]
    node_names=re.findall('(N\d{1,3})',epics_newick)
    for word in node_names:
        epics_newick=epics_newick.replace(word,"",1)
    return epics_newick
 
 
          
#function transforms [0, 0] into [0/0] 
#INPUT - the event
#OUTPUT - the transformed event
def format_event(event):
    formated_event=str(event)
    formated_event=formated_event.replace(', ','/')
    return formated_event



#function writes in file the epicsNewick tree
#INPUT - the path and the tree
#OUTPUT - the file
def write_epicsNewick(epics_newick,name,path):
    f=open(path+name,'w')
    f.write(epics_newick)
    f.close
       

def get_aalist(code):
    #polar
    P='WYTCSNQDEHKR'
    #charged
    c='RKHDE'
    #plus
    p='HKR'
    #minus
    m='DE'
    #hydro
    h='ILVFWYHKTMCGA'
    #aliphtic
    a='ILV'
    #aromatic
    A='FWYH'
    #small
    s='PNTCVAGCSD'
    #tiny
    t='AGCS'
    aalist=''
    for i in code :
        aalist+=eval(i)
    return aalist

       

def add_root(tree):
    root='(R,'
    tree=tree[0:-1]
    newtree=root+tree+');'
    return newtree    
  

       
############MAIN############
#ARG1 : FASTMLpath , ARG2 : Results folder , ARG3 : AAlists used for events (don't put anything
#  for default all AAs) otherwise use combination of letters.     
       
dataname=os.path.split(sys.argv[1][:-1])[1].split('_')[0]
print dataname
FASTMLpath = sys.argv[1]
RESpath= sys.argv[2]
i=0
AAlist=None
if len(sys.argv)==4 :
    i=1
    AAlist=get_aalist(sys.argv[3])
    
   


seqs=read_seqs(FASTMLpath+'seq.marginal.txt') 
ancestortree=read_tree_ancestorformat(FASTMLpath+'tree.ancestor.txt')
events=find_events(seqs,ancestortree,i,AAlist)
epics_newick=make_epicsNewick(FASTMLpath+'tree.newick.txt',events)
rooted_epics=add_root(epics_newick)
write_epicsNewick(rooted_epics,dataname+'_epicsTree',RESpath)
