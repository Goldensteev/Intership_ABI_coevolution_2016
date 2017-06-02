# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 15:11:13 2016

@author: Steven Fletcher
"""

import copy
from Bio import SeqIO
from Bio.Align import AlignInfo
from Bio import pairwise2
import os
import re


#function exctracts seqs of length under fixed threshold form fasta file'
#INPUT - fixed threshold, fastafile and path
#OUTPUT - a file with short seqs, and the len_seq_filtered one  
def Seq_len_filtering(threshold,path,fastafile):
    short_sequences = [] 
    rest_sequences = []
    for record in SeqIO.parse(open(path+fastafile, 'rU'), "fasta") :
        if len(record.seq) < threshold :
            short_sequences.append(record)
        else :
            rest_sequences.append(record)
    print "Found %i short sequences" % len(short_sequences)
    print "Found %i sequences over threshold " % len(rest_sequences)
    short_output_handle = open(path+'short_'+fastafile, "w")
    woshort_output_handle = open(path+'lsf_'+fastafile, "w")
    SeqIO.write(short_sequences, short_output_handle, "fasta")
    SeqIO.write(rest_sequences, woshort_output_handle, "fasta")
    short_output_handle.close()
    woshort_output_handle.close()




#function replace sequence name in MSA file by '>S+#'
#INPUT - the MSA file and the path
#OUTPUT - the renamed sequences MSAfile and the key file     
def change_name_seq(MSAfile,path):
    f=open(path+MSAfile,'r')
    lines=f.readlines()
    f.close()
    f2=open(path+'renamedseqs_'+MSAfile,'w')
    f3=open(path+'seqkeys_'+MSAfile,'w')
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



#function exctracts seqs that are not redundant from redundant datafile (CATH domain list -> CATH nonredundant S100)
#INPUT - non redundant datafile, redundant datafile and datapth
#OUTPUT - a file containing non redundant sequences (all other files remain intact)   
def Seq_redundant_filtering(nonredun_data_file,redun_data_file,datapath) :
    allnonredun_seqs=[]
    nonredun_seqs=[]
    for nonre_record in SeqIO.parse(open(nonredun_data_file, 'rU'), "fasta") :
        allnonredun_seqs.append(nonre_record.id)
    for re_record in SeqIO.parse(open(datapath+redun_data_file, 'rU'), "fasta") :
        if re_record.id in allnonredun_seqs :
            nonredun_seqs.append(re_record)
    nonre_output_handle=open(datapath+'nonred_'+redun_data_file,'w')
    SeqIO.write(nonredun_seqs, nonre_output_handle, "fasta")
    nonre_output_handle.close()
    
    
    
#function merges 2 fasta files and checks for redundancies in them
#INPUT - both fasta file and the datapath
#OUTPUT - A merged fasta file with no redundancies
def Merge_fasta(fastafile1,fastafile2,datapath):
    newfastaid=[]
    newfastaseq=[]    
    for f1record in SeqIO.parse(open(fastafile1,'rU'),'fasta'):
        newfastaid.append(f1record.id)
        newfastaseq.append(f1record)
    for f2record in SeqIO.parse(open(fastafile2,'rU'),'fasta'):
        if f2record.id not in newfastaid :
            newfastaseq.append(f2record)
    mergedf_output_handle=open(datapath+'merged_fasta.txt','w')
    SeqIO.write(newfastaseq, mergedf_output_handle, "fasta")
    mergedf_output_handle.close()



#function checks if the sequences in fasta file are in a database(CATH in this case)
#INPUT - database file and the fasta file
#OUTPUT - A file with all the sequences found by psiblast
def Check_seq_in_database(fastafile,databasefile):
    db={}
    not_inDB=[]
    inDB=[]
    f=open(databasefile)
    lines=f.readlines()
    f.close()
    for line in lines :
        line=line.split()[0][0:5].lower()
        db[line]=None
    for frecord in SeqIO.parse(open(fastafile,'rU'),'fasta') :
        record=(frecord.id.split('_')[0]+frecord.id.split('_')[1]).lower()
        if record in db.keys() :
            inDB.append(frecord)
        else :
            not_inDB.append(frecord)
    return inDB,not_inDB
    
    
#function parse the non redundant pdbdatabase
#INPUT - non redundant pdb file
#OUTPUT - a list of the lines of the nrpdb database
#HOW TO READ LIST
#---------------------------------------------------------------------
# Non-redundant PDB chain set
#---------------------------------------------------------------------
#
# 1: PDB code
# 2: Chain ID
# 3: MMDB ID
#
# 4: Group ID (BLAST pvalue 10e-7)
# 5: Rank (BLAST pvalue 10e-7)
# 6: Representative (1) or not (0) (BLAST pvalue 10e-7)
#
# 7: Group ID (BLAST pvalue 10e-40)
# 8: Rank (BLAST pvalue 10e-40)
# 9: Representative (1) or not (0) (BLAST pvalue 10e-40)
#
# A: Group ID (BLAST pvalue 10e-80)
# B: Rank (BLAST pvalue 10e-80)
# C: Representative (1) or not (0) (BLAST pvalue 10e-40)
#
# D: Group ID (Non-identical sequences)
# E: Rank (Non-identical sequence)
# F: Representative (1) or not (0) (Non-identical sequence)
#
# G: Percentage of Unknown residues
# H: Percentage of Incomplete residues
# I: Percentage of Missing residues
# J: Percentage of Incomplete side-chain residues
# K: Resolution (0.0 if NMR)
# L: No. of chains (subunits) in the PDB entry
# M: No. of heterogens in the PDB entry
# N: No. of different heterogen types in the PDB entry
# O: No. of residues in the chain
# P: Method of coordinate determination (X: X-ray, N: NMR, M: theoretical)
# Q: Acceptable (a) in structural quality or not (n)
#
#---------------------------------------------------------------------------------------------------------------------------
# 1  2       3     4    5 6     7    8 9     A    B C     D    E F    G      H      I      J     K     L   M   N     O P Q
#---------------------------------------------------------------------------------------------------------------------------
def read_nrpdb(nrpdbfile):
    f=open(nrpdbfile,'r')
    nrpdb=f.readlines()
    f.close
    for i in range(len(nrpdb)):
        nrpdb[i]=nrpdb[i].split()
    return nrpdb


#function cleans the pdbseqres (redundant) and makes it non redundant
#INPUT - non redundant pdb file, and the pdb file 
#OUTPUT - a fastafile of the none redundant PDB
def clean_pdbseqres(nrpdb,pdbseqres):
    dic_nrpdb={}
    new_nrpdb=[]
    for record in SeqIO.parse(open(pdbseqres,'rU'),'fasta'):
        dic_nrpdb[record.id.lower()]=record
    for i in nrpdb :
        if i[0].lower()+'_'+i[1].lower() in dic_nrpdb.keys():
            new_nrpdb.append(dic_nrpdb[i[0].lower()+'_'+i[1].lower()])
    output_handle=open('/Users/steven/Desktop/Stage2016/DATA/PDB/NoCATH_lsf_nrpdb_seqres.txt','w')
    SeqIO.write(new_nrpdb, output_handle, "fasta")
    output_handle.close()
    

##NEED TO CLEAN##

#function cleans the redundant CATH file and makes in none redundant using nrpdb
#INPUT - non redundant pdb file an the CATH redundantfile
#OUTPUT - a fastafile of the none redundant PDB
def clean_CATH(nrpdb_seqres,CATHfile):
    dic_nrpdb={}
    cath_in_pdb=[]
    cath_not_pdb=[]
    for record in SeqIO.parse(open(nrpdb_seqres,'rU'),'fasta'):
        dic_nrpdb[(record.id.split('_')[0]+record.id.split('_')[1]).lower()]=record
    for CATHrecord in SeqIO.parse(open(CATHfile,'rU'),'fasta'):
        if CATHrecord.id.split('|')[2].split('/')[0][0:5].lower() in dic_nrpdb.keys():
            cath_in_pdb.append(CATHrecord)                    
        else:
            cath_not_pdb.append(CATHrecord)
#    output_handle=open('/Users/steven/Desktop/Stage2016/DATA/Cytochrome_P450/nredun_S100_CYp450_1.10.630.10.ATOM.fasta.txt','w')
#    SeqIO.write(cath_in_pdb, output_handle, "fasta")
#    output_handle.close()
    return(cath_in_pdb,cath_not_pdb) 



#function cleans the not in cath sequences to get rid of the ones that are not representative compared to the ones in CATH
#INPUT - the CATH sequences in the pdb, those not in and the nrpdb 
#OUTPUT - a list of CATH sequences not in the pdb and without the pdb duplicates
def Filter_CATH_representPDB(cath_in_pdb,cath_not_pdb):
    l=copy.deepcopy(cath_not_pdb)
    for i in range(len(cath_in_pdb)):
        for j in range(len(cath_not_pdb)):
            if cath_in_pdb[i].id.split('|')[2].split('/')[0][0:4].lower() == cath_not_pdb[j].id.split('|')[2].split('/')[0][0:4].lower():
                l[j]=None
    return l


#function cleans the not in cath sequences to get rid of the ones that are not in nrpdb
#INPUT - the CATH sequences in the pdb, the nrpdb 
#OUTPUT - a list of CATH sequences in the pdb but that are not the representative
def Filter_CATH_not_in_representPDB(cath_not_pdb,nrpdb):
    dic_nrpdb={}
    l=copy.deepcopy(cath_not_pdb)
    final_cath_in_pdb=[]
    for i in nrpdb :
        dic_nrpdb[i[0]]=i[1:]
    for j in range(len(cath_not_pdb)) :
        if cath_not_pdb[j]!=None :
            if cath_not_pdb[j].id.split('|')[2].split('/')[0][0:4].upper() in dic_nrpdb.keys():
                 final_cath_in_pdb.append(cath_not_pdb[j])
            else :
                l[j]=None
    return final_cath_in_pdb

 
                                   
#function Compare the CATH domain and the CATHlist.S100 (no 100% id)
#INPUT - CATH domain that we are intersted in and the List of CATH sequences in the domain we are interested in
#OUTPUT - gives in output the list of S100 CATH sequences of the Domain we are interested in 
def get_S100_CATH(CATHdomain,CATHlistS100):
    f=open(CATHlistS100,'r')
    lines=f.readlines()
    f.close()
    lines_dic={}
    for line in lines :
        line=line.split()
        lines_dic[line[0]]=' '+line[1]+'.'+line[2]+'.'+line[3]+'.'+line[4]
    CATHdomainS100=[]
    for record in SeqIO.parse(open(CATHdomain,'rU'),'fasta'):
        if record.id.split('|')[2].split('/')[0] in lines_dic.keys():
            CATHdomainS100.append(record)
    return CATHdomainS100
    


def write_Fasta(path,fasta_filename,Data):
    f=open(path+fasta_filename,'w')
    SeqIO.write(Data,f, "fasta")
    f.close()
    


    
def do_nrpdb_noCATH(nrpdb_representatif,CATHlist):
    f=open(nrpdb_representatif,'r')
    nrpdblines=f.readlines()
    f.close()
    nrdict={}
    for i in nrpdblines :
        j=i.split()[0]+i.split()[1]
        j=j.lower()
        nrdict[j]=i
    f=open(CATHlist,'r')
    cathlist=f.readlines()
    f.close()
    for cathseq in cathlist :
        cathseq=cathseq.split()[0][0:5].lower()
        if cathseq in nrdict.keys() :
            nrdict.pop(cathseq)
    return nrdict
 

   
def write_nrpdb_NoCATH(nrdict):
    f=open('/Users/steven/Desktop/Stage2016/DATA/PDB/nrpdb_represntatifs_noCATH.csv','w')
    for i in nrdict.items() :
        f.write(i[1])
    f.close()



def do_Benchfile(Datafile,lenmax):
    i=0
    f=[]
    for record in SeqIO.parse(open(Datafile,'rU'),'fasta'):
        if i < lenmax :
            f.append(record)
        else :
            break
        i+=1
    write_Fasta('/Users/steven/Desktop/Stage2016/FASTMLbenchmarking/','MSA_'+str(i)+'.fasta.txt',f)



def replace_seq(fastafile) :
    f=open(fastafile,'r')
    lines=f.readlines()
    f.close()
    f=open('/Users/steven/Desktop/Stage2016/DATA/Cytochrome_P450/MSAseq.txt','w')
    for i in lines :
        if i[0]=='>' :
            i=i.replace('-','_')
            i=i.replace('/','_')
            i=i.replace(' ','_')
            i=i.replace('|','_')
            i=i.replace(':','_')
        f.write(i)
    f.close()



def replace_ancestor(ancestorfile):
    f=open(ancestorfile,'r')
    lines=f.readlines()
    f.close()
    f=open('/Users/steven/Desktop/Stage2016/FASTML_runs/CytochromeP450/S100/ances.txt','w')
    for i in lines:
        i=i.replace('-','_')
        f.write(i)
    f.close()

def read_fasta(fastafile):
    seqs=[]
    for record in SeqIO.parse(open(fastafile,'rU'),'fasta'):
        seqs.append(record)
    return seqs



#listdirect=os.listdir('/Users/steven/Desktop/Stage2016/FASTML_runs/ALL_CATH_FAM/')
#listfiles=os.listdir('/Users/steven/Desktop/Stage2016/DATA/CATH_DATA/SequenceBySuperfamily/COMBS_CONSENSUS_TRIM/NR/MSA_MA/')
#renamed=[]
#notin=[]
#for i in listfiles :
#    if 'renamedseqs_' in i :
#        renamed.append(i)    
#f=open('/Users/steven/Desktop/Stage2016/notin.txt')
#l=f.readlines()
#f.close()


for i in l:
    print l
    os.system('scp -r /Users/steven/Desktop/Stage2016/DATA/CATH_DATA/SequenceBySuperfamily/COMBS_CONSENSUS_TRIM/NR/MSA_MA/'+str(l)+' mcarpent@gnome.dsi.upmc.fr:./')


#for i in renamed :
#    if i+'_outDir' not in listdirect :
#        notin.append(i)

#f=open('/Users/steven/Desktop/Stage2016/notin.txt','w')
#for i in notin :        
#    f.write(i+'\n')
#f.close()



#l=os.listdir('/Users/steven/Desktop/Stage2016/FASTML_runs/ALL_CATH_FAM_2/')[1:]
#print len(l)
#for i in l :
#    os.system('mv /Users/steven/Desktop/Stage2016/FASTML_runs/ALL_CATH_FAM/'+i+'/epicsTree.txt '+'/Users/steven/Desktop/Stage2016/FASTML_runs/ALL_CATH_FAM/'+i+'/epicsTree.txt'  )    



#f=open('/Users/steven/Desktop/Stage2016/listFASTML','r')
#listALI=f.readlines()
#f.close()
#print listALI[0]
     
#summary_align = AlignInfo.SummaryInfo('/Users/steven/Desktop/Stage2016/MSA_MA_S100.1.10.8.10_nC.lsf.nrpdb.rep.fasta.txt')
#consensus = summary_align.dumb_consensus()

############MAIN############     
#replace_seq('/Users/steven/Desktop/Stage2016/DATA/Cytochrome_P450/MSA_CYp450_S100_nC.lsf.nrpdb.rep.fasta.txt')  
#replace_ancestor('/Users/steven/Desktop/Stage2016/FASTML_runs/CytochromeP450/S100/tree.ancestor.txt')      
#l=[50,100,150]
#for i in l :
#    do_Benchfile('/Users/steven/Desktop/Stage2016/DATA/Cytochrome_P450/MSA_CYp450_S100_end.fasta.txt',i)    
#nrdict=do_nrpdb_noCATH('/Users/steven/Desktop/Stage2016/DATA/PDB/nrpdb_represntatifs.csv','/Users/steven/Desktop/Stage2016/DATA/CATH_DATA/CathDomainList.txt')    
#write_nrpdb_NoCATH(nrdict)    
#S100=get_S100_CATH('/Users/steven/Desktop/Stage2016/DATA/Cytochrome_P450/CYp450_ALL_1.10.630.10.ATOM.fasta.txt','/Users/steven/Desktop/Stage2016/DATA/CATH_DATA/CathDomainList.S100.txt') 
#write_Fasta('/Users/steven/Desktop/Stage2016/DATA/Cytochrome_P450/','S100_CYp450_1.10.630.10.ATOM.fasta.txt',S100)               
#nrpdb=read_nrpdb('/Users/steven/Desktop/Stage2016/DATA/PDB/nrpdb_represntatifs.csv')
#cath_pdb=clean_CATH('/Users/steven/Desktop/Stage2016/DATA/PDB/lsf_nrpdb_seqres.txt','/Users/steven/Desktop/Stage2016/DATA/Cytochrome_P450/CYp450_CATH_1.10.630.10.ATOM.fasta.txt')        
#print len(cath_pdb[1])
#Get_good_sequences('/Users/steven/Desktop/Stage2016/DATA/Cytochrome_P450/filtered_lenCATHseq.txt','/Users/steven/Desktop/Stage2016/DATA/Cytochrome_P450/CYp450_S100_psib_on_nrpdb.fasta.txt')
#filter_amino_acid('/Users/steven/Desktop/Stage2016/DATA/Cytochrome_P450/','/Users/steven/Desktop/Stage2016/DATA/Cytochrome_P450/CYp450_S100_Lsf_psib_on_nrpdb.fasta.txt',['U'])                
#cath_not_pdb=Filter_CATH_representPDB(cath_pdb[0],cath_pdb[1])
#print Filter_CATH_not_in_representPDB(cath_not_pdb,nrpdb)
#clean_pdbseqres(nrpdb,'/Users/steven/Desktop/Stage2016/DATA/PDB/lsf_pdb_seqres.txt')
#Write_len_Rfile('/Users/steven/Desktop/Stage2016/DATA/Cytochrome_P450/','/Users/steven/Desktop/Stage2016/DATA/Cytochrome_P450/CYp450_S100_psib_on_nrpdb.fasta.txt')          
#change_name_seq('MSA_3.30.1360.30-2226-Aspartate-tRNAligase-mitochondrial.txt','/Users/steven/Desktop/Stage2016/Gyrase_data/')    
#Seq_len_filtering(30,'/Users/steven/Desktop/Stage2016/DATA/PDB/','nrpdb_seqres.txt')
#Seq_redundant_filtering('/Users/steven/Desktop/Stage2016/DATA/CATH_DATA/Cath.DataSet.NonRedundant.S40_overlap_60.v4_0_0.fa.txt','CYp450_1.10.630.10.ATOM.txt','/Users/steven/Desktop/Stage2016/DATA/Cytochrome_P450/')
#Get_psiblast_results('/Users/steven/Desktop/Stage2016/DATA/Cytochrome_P450/psib_nrpdb_S100_CYp450_1.10.630.10.ATOM.txt','/Users/steven/Desktop/Stage2016/DATA/PDB/lsf_nrpdb_seqres.txt','/Users/steven/Desktop/Stage2016/DATA/Cytochrome_P450/')
#Merge_fasta('/Users/steven/Desktop/Stage2016/DATA/Cytochrome_P450/CYp450_psib_on_nrpdb.fasta.txt','/Users/steven/Desktop/Stage2016/DATA/Cytochrome_P450/CYp450_Uc_consensus_psib_on_nrpdb.fasta.txt','/Users/steven/Desktop/Stage2016/DATA/Cytochrome_P450/')
#Uc_CATH=Check_seq_in_database('/Users/steven/Desktop/Stage2016/DATA/Cytochrome_P450/CYp450_goodUc_consensus_psib_on_pdb.fasta.txt','/Users/steven/Desktop/Stage2016/DATA/CATH_DATA/CathDomainList.S100.txt')


    