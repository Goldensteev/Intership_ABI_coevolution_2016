# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 16:10:37 2016

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
    

def main1():
    listfams=glob.glob('/Users/steven/Desktop/Stage2016/SpecificFams/new/MSAcleaned/MSA_*')
    badcolumns={}
    for i in listfams :
        fam,seqs=parseMSAfile(i)
        badcolumns[fam]=[]
        poslist=[]
        for pos in range(len(seqs[0].seq)+1):
            poslist.append(pos+1)
        for seq in seqs :
            for posaa in range(1,len(seq.seq)+1) :
                if seq.seq[posaa-1]=='-' and posaa in poslist:
                    poslist.remove(posaa)
                    print 'removed '+str(posaa)
                
        badcolumns[fam]=poslist
    print len(badcolumns['3.10.200.10'])

def MAIN():
    f=open('/Users/steven/Desktop/Stage2016/SpecificFams/new/ExempleFamilles/OK/2.40.180.10/AllRes.txt','r')
    lines=f.readlines()
    print lines[1]
    f.close()

MAIN()               
            