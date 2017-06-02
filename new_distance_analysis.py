# -*- coding: utf-8 -*-
"""
Created on Tue May 24 15:45:33 2016

@author: steven
"""

import ATOMStoCOMBS as AtC
import os
from Bio import SeqIO
import operator


def do_distmatrix(listallCATH):
    for i in listallCATH :
        print i
        os.system('/Users/steven/Desktop/Stage2016/CalcDistCA/CalcDistCA/./computeCADist -i /Users/steven/Desktop/Stage2016/DATA/CATH_DATA/dompdb/'+i+' > /Users/steven/Desktop/Stage2016/dist_struct_Mat/'+i)
 

       
def parseMSAfile(MSApath,MSAfile):
    seqs=[]
    fam=MSAfile.split('_')[2]
    for record in SeqIO.parse(open(MSApath+MSAfile, 'rU'), "fasta") :           
        seqs.append(record)
    return fam,seqs



def get_distmat(seqname_distmat):
    test=TESTexist('/Users/steven/Desktop/Stage2016/dist_struct_Mat/'+seqname_distmat)
    if test==False:
        return False
    f=open('/Users/steven/Desktop/Stage2016/dist_struct_Mat/'+seqname_distmat)
    distmat=f.read().splitlines()
    for i in range(len(distmat)) :
        distmat[i]=distmat[i].split()
    f.close()
    distmat.pop(0)
    distmat.pop(-1)
    for i in distmat :
        i.pop(0)
    return distmat
    
    
    
def get_epicsmat(famname,epicsTYPE):
    test=TESTexist('/Users/steven/Desktop/Stage2016/epicsMAT/'+epicsTYPE+'/'+famname+'_mat'+epicsTYPE)
    if test==False:
        return False
    f=open('/Users/steven/Desktop/Stage2016/epicsMAT/'+epicsTYPE+'/'+famname+'_mat'+epicsTYPE)
    epicsmat=f.read().splitlines()
    for i in range(len(epicsmat)) :
        epicsmat[i]=epicsmat[i].split()
    f.close()
    return epicsmat
    
    
    
def get_epicsmatpval(famname,epicsTYPE):
    test=TESTexist('/Users/steven/Desktop/Stage2016/epicsMAT/'+epicsTYPE+'/'+famname+'_matPVAL')
    if test==False:
        return False    
    f=open('/Users/steven/Desktop/Stage2016/epicsMAT/'+epicsTYPE+'/'+famname+'_matPVAL')
    epicsmat=f.read().splitlines()
    for i in range(len(epicsmat)) :
        epicsmat[i]=epicsmat[i].split()
    f.close()
    return epicsmat    



def parseSTRIDEfiles(STRIDEfile_list):
    STRIDE_dict={}
    stride_list=[]
    f=open(STRIDEfile_list,'r')
    stride_list=f.read().splitlines()
    f.close()
    for i in stride_list :
        print i
        f=open('/Users/steven/Desktop/Stage2016/DATA/CATH_DATA/Stride_dompdb_out/'+i,'r')
        temp=f.readlines()
        for j in range(len(temp)) :
            temp[j]=temp[j].split()    
        STRIDE_dict[i]=temp
    f.close()
    return STRIDE_dict    



def get_2nd_struct(seqname,posAA1,posAA2):
    f=open('/Users/steven/Desktop/Stage2016/DATA/CATH_DATA/Stride_dompdb_out/'+seqname,'r')
    temp=f.readlines()
    for j in range(len(temp)) :
        temp[j]=temp[j].split()
    structAA1=0
    structAA2=0
    for i in temp:
        if i[0]=='ASG' :
            if int(i[4])==posAA1 :
                structAA1=i[6]
            if int(i[4])==posAA2 :
                structAA2=i[6]
    return structAA1,structAA2



def get_2nd_struct_bis(seqname,posAA):
    f=open('/Users/steven/Desktop/Stage2016/DATA/CATH_DATA/Stride_dompdb_out/'+seqname,'r')
    temp=f.readlines()
    for j in range(len(temp)) :
        temp[j]=temp[j].split()
    structAA=0
    enfouiAA=0
    for i in temp:
        if i[0]=='ASG' :
            if int(i[4])==posAA :
                structAA=i[6]
                enfouiAA=i[9]
    return structAA,enfouiAA



def getATOMS():
    ATOMSseqs={}
    for record in SeqIO.parse(open('/Users/steven/Desktop/Stage2016/ATOMS.aligned.txt','rU'), "fasta") :
        ATOMSseqs[record.id.split('|')[2].split('/')[0]]=record
    return ATOMSseqs
 
       
        
def getnbGAPS(seq,posAA):
    nbgaps=0
    for pos in range(len(seq)) :
        if seq[pos]=='-' :
            nbgaps+=1
        if pos >= posAA :
            return nbgaps


def writeRES(family,familyname):
    #f=open('/Volumes/DDEXT/RES2/ '+familyname,'w')
    f=open('/Users/steven/Desktop/Stage2016/RES/ '+familyname,'w')
    f.write('SEQ\t  option\t  AA1\t  posAA1\t  struAA1\t  enfouiAA1\t  AA2\t  posAA2\t  struAA2\t  enfouiAA1\t  distpair\t  nbepics\t  pVAL\n')
    for sequence in family.keys() :
        for option in family[sequence].keys() : 
            for res in family[sequence][option]:
                f.write(sequence+'\t '+option+'\t ')
                for i in res :
                    f.write(str(i)+'\t ')
                f.write('\n')
            f.write('\n')
    f.close()      
    



def TESTexist(name):
    test=os.path.isfile(name)
    return test