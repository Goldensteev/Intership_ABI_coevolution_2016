# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 16:44:08 2016

@author: steven
"""

from Bio import SeqIO
from numpy import std

#function does the standard deviation according the lengths of the sequences
#INPUT - fastafile 
#OUTPUT - the standard deviation
def get_standard_deviation(CATHfile):
    lenseqs=[]
    for record in SeqIO.parse(open(CATHfile,'rU'),'fasta'):
        lenseqs.append(len(record))
    return std(lenseqs)  

#function reads a fastafile and returns the sequences and all other info
#INPUT - fastafile 
#OUTPUT - the sequences of the fastafile
def read_fasta(fastafile):
    CATHseqs=[]
    for record in SeqIO.parse(open(fastafile,'rU'),'fasta'):
        if 'cath' in record.id :
            CATHseqs.append(record)
        else :
            BIGseq=record
    return CATHseqs,BIGseq


def trim_lenghty_sequences(seqs,consesus):
    cindex=[]
    count=0
    for pos in consesus.seq :
        if pos !='-':
            cindex.append(count)
        count+=1   
    cindex=(cindex[0],cindex[-1])
    for s in seqs :
        s.seq=s.seq[cindex[0]:cindex[-1]]
    return s

def trim_lenghty_seq(CATHseqs,BIGseq):
    aapos=[]
    for seq in CATHseqs :
        for pos in range(len(seq.seq)) :
            if seq.seq[pos] !='-' and pos not in aapos:
                aapos.append(pos)
    print sorted(aapos)                       
    museq=BIGseq.seq.tomutable() 
    count=0           
    for pos in range(len(BIGseq.seq)):
        if pos not in aapos : 
            museq.pop(pos-count)
            count+=1
        elif museq[pos-count]=='-':
            museq.pop(pos-count)
            count+=1
    print museq
    BIGseq.seq=museq.toseq()            
    return BIGseq


print get_standard_deviation('/Users/steven/Desktop/Stage2016/DATA/CATH_DATA/SequenceBySuperfamily/COMBS_TEST/NR/nonred_1.10.8.10.COMBS.v4.0.0.fasta.txt')
s,B=read_fasta('/Users/steven/Desktop/Stage2016/outMUprofile.txt')
#print trim_lenghty_seq(s,B)