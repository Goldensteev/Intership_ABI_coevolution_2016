# -*- coding: utf-8 -*-
"""
Created on Fri May 27 16:14:00 2016

@author: steven
"""

import re

def read_description_file(description_file):
    f=open(description_file,'r')
    description=f.read().splitlines()
    f.close()
    return description



def parse_description_file(description):
    des={}
    for i in range(len(description)) :
        if 'DOMAIN' in description[i] :
            j=description[i].split()
            des[j[1]]=''
        if 'SOURCE' in description[i]:
            des[j[1]]=des[j[1]]+description[i][10:]
            #des[j[1]]=des[j[1]].split('.')[0]
    return des



def get_seqkeys(SeqkeysFile):
    f=open(SeqkeysFile,'r')
    lines=f.read().splitlines()
    f.close()
    seqkeys={}
    for i in range(len(lines)):
        lines[i]=lines[i].split('=')
        lines[i][-1]=lines[i][-1][3:]
        if '>cath' in lines[i][0] :
            lines[i][0]=lines[i][0].split('|')[2].split('/')[0]
        seqkeys[lines[i][-1]]=lines[i][0]
    return seqkeys



def replace_name(tree,des,seqkeys):
    f=open(tree,'r')
    newnamestree=f.read().splitlines()
    f.close()
    node_names=re.findall('(S\d{1,3})',newnamestree[0])
    for i in node_names:
        if seqkeys[i] in des:
            j=des[seqkeys[i]]
            print j+'\n'
            newnamestree[0]=newnamestree[0].replace(i,j,1)
        else :
            newnamestree[0]=newnamestree[0].replace(i,seqkeys[i],1)
    return newnamestree[0]



def write_epicsNewick(epics_newick,name,path):
    f=open(path+name,'w')
    f.write(epics_newick)
    f.close



def add_root(tree):
    root='(R,'
    tree=tree[0:-1]
    newtree=root+tree+':0);'
    return newtree    


desc=read_description_file('/Users/steven/Desktop/Stage2016/DATA/CATH_DATA/CathDomainDescriptionFile.v4.0.0')
des=parse_description_file(desc)
seqkeys=get_seqkeys('/Users/steven/Desktop/Stage2016/SpecificFams/new/3.10.200.10/seqkeys_MSA_3.10.200.10.txt')
tree=replace_name('/Users/steven/Desktop/Stage2016/SpecificFams/new/newickTrees/3.10.200.10_tree.newick.txt',des,seqkeys)
write_epicsNewick(tree,'3.10.200.10.test.txt','/Users/steven/Desktop/Stage2016/SpecificFams/new/')
