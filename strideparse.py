# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 12:42:14 2016

@author: steven
"""
import glob
import os
from Bio import SeqIO

def parseMSAfile(MSA):
    seqs=[]
    m=os.path.split(MSA)
    fam=m[1].split('_')[1]
    for record in SeqIO.parse(open(MSA, 'rU'), "fasta") :           
        seqs.append(record)
    return fam,seqs
    

def parsestruct(seqstruct):
    count_list_struct={}
    for fam in seqstruct.keys():
        count_list_struct[fam]={}
        for seq in seqstruct[fam].keys():
            for struct in seqstruct[fam][seq]:
                if struct not in count_list_struct[fam].keys() :
                    count_list_struct[fam][struct]=1
                else :
                    count_list_struct[fam][struct]+=1
    return count_list_struct

def get_2nd_struct(seqname):
    f=open('/Users/steven/Desktop/Stage2016/DATA/CATH_DATA/Stride_dompdb_out/'+seqname,'r')
    temp=f.readlines()
    for j in range(len(temp)) :
        temp[j]=temp[j].split()
    structseq=[]
    for i in temp:
        if i[0]=='ASG' :
            structAA=i[6]
            structseq.append(structAA)        
    return structseq
    

listfams=glob.glob('/Users/steven/Desktop/Stage2016/SpecificFams/new/MSAcleaned/MSA_*')
seqstruct={}
for i in listfams :
    fam,seqs=parseMSAfile(i)
    seqstruct[fam]={}
    for seq in seqs :
        if 'cath' in seq.id :
            for aa in seq.seq:
                seqname=seq.id.split('_')[4]
                seqstruct[fam][seqname]=get_2nd_struct(seqname)
print parsestruct(seqstruct)