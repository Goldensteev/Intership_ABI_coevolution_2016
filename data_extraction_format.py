# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 11:05:55 2016

@author: steven
"""
from Bio import SeqIO
from numpy import std
import os
import sys

#function writes a fastafile with Biopython module
#INPUT - datapath, Filename and the List withe the Data 
#OUTPUT - A file with all the sequences found by psiblast
def write_Fasta(DATApath,fasta_filename,Data):
    f=open(DATApath+fasta_filename,'w')
    SeqIO.write(Data,f, "fasta")
    f.close()


#function reads psiblast file and extracts the list of ids of the hits then puts all the sequences in fastafile
#INPUT - psiblast file, the databasefile and the data path
#OUTPUT - A file with all the sequences found by psiblast
def Get_psiblast_results(PDBpath,DATApath,psibfile,pdbfile,dataname,database):
    psib_ids=[]
    pdb_pblasted_seqs=[]
    for psib_record in SeqIO.parse(open(DATApath+psibfile,'rU'),'fasta'):
        psib_ids.append(psib_record.id)
    for pdb_record in SeqIO.parse(open(PDBpath+pdbfile,'rU'),'fasta'):
        if pdb_record.id in psib_ids:
            pdb_pblasted_seqs.append(pdb_record)
    write_Fasta(DATApath,dataname+'_'+database+'.fasta.txt',pdb_pblasted_seqs)
  
   
   
#function merges 2 fasta files and checks for redundancies in them
#INPUT - both fasta file and the datapath
#OUTPUT - A merged fasta file with no redundancies
def Merge_fasta(fastafile1,fastafile2,DATApath,RESpath,dataname):
    newfastaid=[]
    newfastaseq=[]    
    for f1record in SeqIO.parse(open(DATApath+fastafile1,'rU'),'fasta'):
        newfastaid.append(f1record.id)
        newfastaseq.append(f1record)
    for f2record in SeqIO.parse(open(RESpath+fastafile2,'rU'),'fasta'):
        if f2record.id not in newfastaid :
            newfastaseq.append(f2record)
    write_Fasta(RESpath,dataname+'_merged.fasta.txt',newfastaseq)
 

   
#function writes a file for the Rscript seqID / seqLength
#INPUT - path and the fasta file 
#OUTPUT - writes a file with the format : Column 1 seqID  Column 2 seqLength  
def Write_len_Rfile(DATApath,datafile,dataname):
    f=open(DATApath+'R_all_'+dataname+'_lenCATHseq.txt','w')
    for record in SeqIO.parse(open(DATApath+datafile,'rU'),'fasta'):
        f.write(str(record.id)+' '+str(len(record.seq))+'\n')
    f.close()
    
    
    
#function gets the list of the good length sequences from the R script and writes the Fastafile with the list
#INPUT - The R file that has the list of sequence names and the CATHfile that is in fasta form
#OUTPUT - Fasta file of the seuences of the right length  
def Get_good_sequences_from_Rlenfile(R_good_seqfile,datafile,DATApath,dataname,database):
    f=open(R_good_seqfile,'r')
    lines=f.readlines()
    f.close()
    good_seqs=[]
    for record in SeqIO.parse(open(datafile,'rU'),'fasta'):
        if record.id+'\n' in lines :
            good_seqs.append(record)
    write_Fasta(DATApath,dataname+'_LSF_'+database+'.fasta.txt',good_seqs)        
    return good_seqs



#function filters out the non recognised AAs from a fasta file
#INPUT - Fastafile and the list of Non AAs
#OUTPUT - writes a file fasta file with sequences with the non recognized aa by '?' (Better for FASTML)  
def filter_amino_acid(DATApath,Fastafile,NOTaalist,dataname,database):
    f=open(Fastafile,'r')
    lines=f.readlines()
    f.close()
    filteredfasta=[]
    for line in lines :
        if line[0] != '>' :
            for Naa in NOTaalist :
                line=line.replace(Naa,'?')
        filteredfasta.append(line)        
    f=open(DATApath+dataname+'_FAA.LSF_'+database+'.fasta.txt','w')
    for i in filteredfasta :
        f.write(i)
    f.close()



#function replace sequence name in MSA file by '>S+#'
#INPUT - the MSA file and the path
#OUTPUT - the renamed sequences MSAfile and the key file     
def change_name_seq(MSAfile,MSApath):
    f=open(MSApath+MSAfile,'r')
    lines=f.readlines()
    f.close()
    f2=open(MSApath+'renamedseqs_'+MSAfile,'w')
    f3=open(MSApath+'seqkeys_'+MSAfile,'w')
    i=0
    for line in lines :
        if line[0]=='>':
            i+=1
            f3.write(line[:-1]+'  =  >S'+str(i)+'\n')
            f2.write('>S'+str(i)+'\n')
        else :
            f2.write(line)
    f2.close()
    f3.close()  



#function checks the # of sequences in a file
#INPUT - the Datapath and the fastafile
#OUTPUT - the number of file 
def check_nbseqs(DATApath,fastafile):
    seqs=[]
    for record in SeqIO.parse(open(DATApath+fastafile,'rU'),'fasta'):
        seqs.append(record)
    return len(seqs)



#function does the standard deviation according the lengths of the sequences
#INPUT - fastafile 
#OUTPUT - the standard deviation
def get_standard_deviation(CATHfile):
    lenseqs=[]
    for record in SeqIO.parse(open(CATHfile,'rU'),'fasta'):
        lenseqs.append(len(record))
    return std(lenseqs)  



#function checks the length of the sequences and if they are over the threshold
#INPUT - the path, the fastafile and the threshold
#OUTPUT - NONE      
def check_len_trim(DATApath,fastafile,threshold):
    longseqs=[]
    nbsq=0
    for record in SeqIO.parse(open(DATApath+fastafile,'rU'),'fasta'):
        if len(record.seq) > (threshold):
            longseqs.append(record)
            write_Fasta(DATApath+'TRIM/','S'+str(nbsq)+'.txt',record)    
            nbsq+=1                                 
    return nbsq        



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



#function trims the lengthy sequences according to the alignment with the consensus
#INPUT - the long sequences and the consensus sequences
#OUTPUT - a list of trimmed long sequences 
def trim_lenghty_seq(CATHseqs,BIGseq):
    aapos=[]
    for seq in CATHseqs :
        for pos in range(len(seq.seq)) :
            if seq.seq[pos] !='-' and pos not in aapos:
                aapos.append(pos)                       
    museq=BIGseq.seq.tomutable() 
    count=0           
    for pos in range(len(BIGseq.seq)):
        if pos not in aapos : 
            museq.pop(pos-count)
            count+=1
        elif museq[pos-count]=='-':
            museq.pop(pos-count)
            count+=1
    BIGseq.seq=museq.toseq()            
    return BIGseq



#function runs all the command lines in order to trim
#INPUT - The Result path, the DATApath ,the CATHfile, the psiblast file the avergage length of the CATH seqs
#OUTPUT - NONE
def do_trimming(RESpath,DATApath,datafile,psibfasta,averagelen):
    trimmedseqs=[]
    os.mkdir(RESpath+'TRIM/')
    os.system('/Users/steven/Desktop/Stage2016/MUSCLE/muscle3.8.31_i86darwin64 -in '+DATApath+datafile+' -out '+RESpath+'TRIM/MU_'+datafile)
    std=get_standard_deviation(DATApath+datafile)
    nbseq=check_len_trim(RESpath,psibfasta,std+averagelen)
    for i in range(nbseq):
        os.system('/Users/steven/Desktop/Stage2016/MUSCLE/muscle3.8.31_i86darwin64 -profile -in1 '+RESpath+'TRIM/MU_'+datafile+' -in2 '+RESpath+'TRIM/S'+str(i)+'.txt -out '+RESpath+'TRIM/out'+str(i)+'.txt')
        cathseqs,bigseq=read_fasta(RESpath+'TRIM/out'+str(i)+'.txt')
        trimmedseqs.append(trim_lenghty_seq(cathseqs,bigseq))
    write_Fasta(RESpath,dataname+'_LSF_'+database+'.fasta.txt',trimmedseqs)
    
    
    
#function calculates mean of seq length
#INPUT - The fastafile
#OUTPUT - the mean of the length
def get_average_length(CATHfile):
    mean=0
    c=0
    for record in SeqIO.parse(open(CATHfile,'rU'),'fasta'):
        mean+=len(record.seq)
        c+=1
    mean/=c
    return mean
    
    
    
    
    
####################################MAIN####################################
#ARG1 : DATAFILE
DATA=os.path.split(sys.argv[1])
DATApath=DATA[0]+'/'
datafile=DATA[1]
a=datafile.split('_')[1]
dataname='S100.'+a.split('.')[0]+'.'+a.split('.')[1]+'.'+a.split('.')[2]+'.'+a.split('.')[3]
print dataname

#ARG2 : DATABASE
PDB=os.path.split(sys.argv[2])
PDBpath=PDB[0]+'/'
pdbfile=PDB[1]
database='nC.lsf.nrpdb.rep'

#ARG3
MSA=sys.argv[3]

#AAs not recognized by MSA programs
NOTaalist=['U','O','X','B','Z','J']
aalist=['ARNDCEQGHILKMFPSTWYV']

RESpath=DATApath+dataname+'/'
MSApath='/Users/steven/Desktop/Stage2016/DATA/CATH_DATA/SequenceBySuperfamily/COMBS_CONSENSUS_TRIM/NR/MSA_MA/'
psibfile='psib_'+dataname+'_on_'+database+'.txt'
Rinfile=RESpath+'R_all_'+dataname+'_lenCATHseq.txt'
Routfile=RESpath+'R_'+dataname+'_filtered_lenCATHseq.txt'
os.mkdir(DATApath+dataname)

####CALL PSIBLAST####
os.system('psiblast -query '+DATApath+datafile+' -subject '+PDBpath+pdbfile+' -evalue=0.000001 -out '+RESpath+psibfile)
Get_psiblast_results(PDBpath,RESpath,psibfile,pdbfile,dataname,database)
averagelen=get_average_length(DATApath+datafile)
do_trimming(RESpath,DATApath,datafile,dataname+'_'+database+'.fasta.txt',averagelen)
Merge_fasta(datafile,dataname+'_LSF_'+database+'.fasta.txt',DATApath,RESpath,dataname)
lenfile=check_nbseqs(RESpath,dataname+'_merged.fasta.txt')
if lenfile > 1:
    #Write_len_Rfile(RESpath,dataname+'_Psib_merged.fasta.txt',dataname)
####CALL R SCRIPT####
    #os.system('Rscript --vanilla Filter_len_seq.R '+Rinfile+' '+Routfile)
    #Get_good_sequences_from_Rlenfile(Routfile,RESpath+dataname+'_Psib_merged.fasta.txt',RESpath,dataname,database)
    filter_amino_acid(RESpath,RESpath+dataname+'_merged.fasta.txt',NOTaalist,dataname,database)
####CALL MSA(MAFFT/CLOSTALO/MUSCLE####
    if MSA == 'ma' :
        os.system('mafft '+RESpath+dataname+'_FAA.LSF_'+database+'.fasta.txt > '+MSApath+'MSA_MA_'+dataname+'_'+database+'.fasta.txt')
        change_name_seq('MSA_MA_'+dataname+'_'+database+'.fasta.txt',MSApath)
    else :    
        os.system('/Users/steven/Desktop/Stage2016/MUSCLE/muscle3.8.31_i86darwin64 -in '+RESpath+dataname+'_FAA.LSF_'+database+'.fasta.txt'+' -out '+DATApath+'MSA_MU/'+'MSA_MU_'+dataname+'_'+database+'.fasta.txt')
        change_name_seq('MSA_MU_'+dataname+'_'+database+'.fasta.txt',MSApath)
####FASTML####
############################################################################