# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 12:08:10 2016

@author: steven
"""

from Bio import SeqIO
import glob
import os
from numpy import std



def get_CATH_infos(cathfile):
    infos={}
    lenseqs=[]
    seqs=[]
    for record in SeqIO.parse(open(cathfile, 'rU'), "fasta") :
        seqs.append(record)           
        lenseqs.append(len(record.seq))
    moy=0
    for i in lenseqs :
        moy=moy+i
    moy=moy/len(lenseqs)    
    infos['len_domaines']=moy
    infos['nb_seqs']=len(seqs)
    return infos

def write_info(CATHinfos):
    f=open('/Users/steven/Desktop/Stage2016/SpecificFams/new/CATHinfoMSA.csv','w')
    f.write('CATHfamily'+'\t'+'nb_seqs'+'\t'+'len_domaines'+'\n')
    for i in CATHinfos.keys() :
        f.write(i+'\t'+str(CATHinfos[i]['nb_seqs'])+'\t'+str(CATHinfos[i]['len_domaines'])+'\n')
    f.close()

def MAIN():
    listCATH=glob.glob('/Users/steven/Desktop/Stage2016/DATA/CATH_DATA/SequenceBySuperfamily/COMBS_CONSENSUS_TRIM/NR/MSA_MA/MSA*')
    CATHinfos={}    
    for i in listCATH :
        dataname=os.path.split(i)[1].split('_')[1][0:11]
        CATHinfos[dataname]=get_CATH_infos(i)
    write_info(CATHinfos)

    
MAIN()