# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 10:22:16 2016

@author: steven
"""

from Bio import SeqIO
import os 


#function exctracts seqs of length under fixed threshold form fasta file'
#INPUT - fixed threshold, fastafile and path
#OUTPUT - a file with short seqs, and the len_seq_filtered one  
def Seq_len_filtering(threshold,DATApath,fastafile):
    short_sequences = [] 
    rest_sequences = []
    for record in SeqIO.parse(open(DATApath+fastafile, 'rU'), "fasta") :
        if len(record.seq) < threshold :
            short_sequences.append(record)
        else :
            rest_sequences.append(record)
    DATApath=DATApath+'LSF/'
    write_Fasta(DATApath,'lsf_'+fastafile,rest_sequences)
    
    
    
#function writes a fastafile with Biopython module
#INPUT - datapath, Filename and the List withe the Data 
#OUTPUT - A file with all the sequences found by psiblast
def write_Fasta(DATApath,fasta_filename,Data):
    f=open(DATApath+fasta_filename+'.fasta.txt','w')
    SeqIO.write(Data,f,'fasta')
    f.close()
 


#function parses the CATH S100 fasta list
#INPUT - the CATHS100fastafile
#OUTPUT - a list of the ID sequences found in the file
def readS100file(nonredun_data_file):
    allnonredun_seqs=[]
    for nonre_record in SeqIO.parse(open(nonredun_data_file, 'rU'), "fasta") :
        allnonredun_seqs.append(nonre_record.id)
    return allnonredun_seqs    
 

       
#function exctracts seqs that are not redundant from redundant datafile (CATH domain list -> CATH nonredundant S100)
#INPUT - non redundant datafile, redundant datafile and datapth
#OUTPUT - a file containing non redundant sequences (all other files remain intact)   
def Seq_redundant_filtering(allnonredun_seqs,redun_data_file,datapath) :
    nonredun_seqs=[]
    for record in SeqIO.parse(open(datapath+redun_data_file, 'rU'), "fasta") :
        if record.id in allnonredun_seqs :
            nonredun_seqs.append(record)
    datapath='/Users/steven/Desktop/Stage2016/SpecificFams/new/'        
    write_Fasta(datapath+'NR/','nonred_'+redun_data_file,nonredun_seqs)





#listfiles=os.listdir('/Users/steven/Desktop/Stage2016/DATA/CATH_DATA/SequenceBySuperfamily/COMBS_CONSENSUS_TRIM/NR/')[2:]
CATHS100=readS100file('/Users/steven/Desktop/Stage2016/DATA/CATH_DATA/CathDomainSeqs.S100.COMBS.v4.0.0.txt')
#for fastafile in listfiles:
fastafile='nonred_3.40.50.410.COMBS.v4.0.0.fasta.txt'
Seq_len_filtering(30,'/Users/steven/Desktop/Stage2016/SpecificFams/new/',fastafile)
#Seq_redundant_filtering(CATHS100,fastafile,'/Users/steven/Desktop/Stage2016/SpecificFams/new/' )
os.system('python data_extraction_format.py /Users/steven/Desktop/Stage2016/DATA/CATH_DATA/SequenceBySuperfamily/COMBS_CONSENSUS_TRIM/NR/'+fastafile+' /Users/steven/Desktop/Stage2016/DATA/PDB/NoCATH_lsf_nrpdb_representative_seqres.fasta.txt ma')
       