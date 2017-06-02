# -*- coding: utf-8 -*-
"""
Created on Tue May  3 13:51:56 2016

@author: steven
"""

from Bio import SeqIO
from Bio import pairwise2
import os
import numpy as np


def get_FASTA(FASTAfile) :
    seqs=[]
    for record in SeqIO.parse(open(FASTAfile, 'rU'), "fasta") :           
        seqs.append(record)
    return seqs




def get_FASTAbis(FASTAfile,seqname) :
    for record in SeqIO.parse(open(FASTAfile, 'rU'), "fasta") :
        if record.id == seqname :           
            return record



def do_pairwise_alignment(ATOM,COMBS):
    alignments=[]
    f=open('/Users/steven/Desktop/Stage2016/ATOMS.aligned.txt','w')
    for i in range(len(ATOM)):
        print ATOM[i].id
        #alignments.append(pairwise2.align.globalxx(ATOM[i].seq,COMBS[i].seq))
        f.write('>'+ATOM[i].id+'\n')
        f.write(str(pairwise2.align.globalxx(ATOM[i].seq,COMBS[i].seq)[0][0])+'\n')
    f.close()
    return alignments
 


def get_distmat(seqname_distmat):
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


 
def make_fullmat(a):
    for i in range(len(a)):
        for j in range(i+1):
            a[i].insert(0,0)
    a.insert(len(a),[0]*(len(a)+1))
    return a
    


def make_symetricalmat(a):
    for i in range(len(a)):
        for j in range(len(a)):
            a[j][i]=a[i][j]
    return a    


def reformat_mat_better(a,b):
    listgaps=[]
    for i in range(len(b)):
        if b[i]=='-':
            a.insert(i,[-1]*len(b))
            listgaps.append(i)
            for j in range(len(a)):
                if j not in listgaps :
                    a[j].insert(i,-1)
    return a

  
  
def reformat_mat(oldmat,seq):
    nbgaps=0
    oldmat=np.asarray(oldmat)
    newmat=np.zeros((len(seq),len(seq)))
    for i in range(len(seq)):
        if seq[i]=='-':
            newmat[i]=-1
            newmat[:,i]=-1
            nbgaps+=1
        else :
            nbgapsbis=0
            for j in range(len(seq)) :
                if seq[j]=='-':
                    newmat[i][j]=-1
                    nbgapsbis+=1
                else :
                    newmat[i][j]=oldmat[i-nbgaps][j-nbgapsbis]
        newmat[i][i]=0
    return newmat


   
def write_Mat_file(reformed_mat,seqname):
    f=open('/Users/steven/Desktop/Stage2016/dist_struct_mat_ATOM_new/'+seqname,'w')
    for i in range(len(reformed_mat)):
        if i != 0 :
            f.write('\n')
        for j in reformed_mat[i]:
            f.write(str(j)+'\t ')
    f.close()
  
  

def parse_CATHlist(CATHlist):
    cathlist=[]
    cathdict={}
    f=open(CATHlist,'r')
    cathlist=f.readlines()
    for i in range(len(cathlist)) :
        cathlist[i]=cathlist[i].split()
        famname=cathlist[i][1]+'.'+cathlist[i][2]+'.'+cathlist[i][3]+'.'+cathlist[i][4]
        if famname not in cathdict.keys() :
            cathdict[famname]=[]
            cathdict[famname].append(cathlist[i][0])
        else :
            cathdict[famname].append(cathlist[i][0])
    f.close()
    return cathdict



def do_mafft_addfull_keeplength(seq,cathdict,list_distmat):
    if len(seq.id.split('|'))!=1 :
        seq.id=seq.id.split('|')[2].split('/')[0]
    seq.description=seq.id
    f=open('/Users/steven/Desktop/Stage2016/temp','w')
    SeqIO.write(seq, f, "fasta")
    f.close()
    alignedseq=None
    for key, value in cathdict.iteritems():
        for seqname in value :
            if seqname == seq.id :
                famname=key
    if 'MSA_MA_S100.'+famname+'_nC.lsf.nrpdb.rep.fasta.txt' in os.listdir('/Users/steven/Desktop/Stage2016/DATA/CATH_DATA/SequenceBySuperfamily/COMBS_CONSENSUS_TRIM/NR/MSA_MA/'):
        os.system('mafft --keeplength --addfull /Users/steven/Desktop/Stage2016/temp /Users/steven/Desktop/Stage2016/DATA/CATH_DATA/SequenceBySuperfamily/COMBS_CONSENSUS_TRIM/NR/MSA_MA/MSA_MA_S100.'+famname+'_nC.lsf.nrpdb.rep.fasta.txt > /Users/steven/Desktop/Stage2016/temp_out'+seq.id)
        alignedseq=get_FASTAbis('/Users/steven/Desktop/Stage2016/temp_out'+seq.id,seq.id)
        os.remove('/Users/steven/Desktop/Stage2016/temp_out'+seq.id)
    os.remove('/Users/steven/Desktop/Stage2016/temp')
    return alignedseq
    

      

def MAIN():
    list_distmat=os.listdir('/Users/steven/Desktop/Stage2016/dist_struct_Mat/')[1:]
    cathdict=parse_CATHlist('/Users/steven/Desktop/Stage2016/DATA/CATH_DATA/CathDomainList.S100.txt')
    for record in SeqIO.parse(open('/Users/steven/Desktop/Stage2016/ATOMS.aligned.txt','rU'), "fasta") :
        seq=record
        seqname=seq.id.split('|')[2].split('/')[0]
        alignedseq=do_mafft_addfull_keeplength(seq,cathdict,list_distmat)
        if alignedseq == None :
            continue
        else :
            print seqname
            distmat=get_distmat(seqname)
            print 'get_distmat DONE'
            distmat=make_fullmat(distmat)
            print 'make_fullmat DONE'
            distmat=make_symetricalmat(distmat)
            print 'make_symetricalmat DONE'
            distmat=reformat_mat(distmat,alignedseq.seq)
            print 'reformat_mat DONE'
            write_Mat_file(distmat,seqname)
            

    
        
    
#MAIN()

#ATOM=get_FASTA('/Users/steven/Desktop/Stage2016/DATA/CATH_DATA/CathDomainSeqs.S100.ATOM.v4.0.0.txt')
#COMBS=get_FASTA('/Users/steven/Desktop/Stage2016/DATA/CATH_DATA/CathDomainSeqs.S100.COMBS.v4.0.0.txt')    
#do_pairwise_alignment(ATOM,COMBS)