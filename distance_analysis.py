# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 15:12:59 2016

@author: steven
"""
import ATOMStoCOMBS as AtC
import os
from Bio import SeqIO
import operator
import glob 
import copy


def do_distmatrix(listallCATH):
    for i in listallCATH :
        os.system('/Users/steven/Desktop/Stage2016/CalcDistCA/CalcDistCA/./computeCADist -i /Users/steven/Desktop/Stage2016/DATA/CATH_DATA/dompdb/'+i+' > /Users/steven/Desktop/Stage2016/dist_struct_Mat/'+i)
 

       
def parseMSAfile(MSA):
    seqs=[]
    m=os.path.split(MSA)
    fam=m[1].split('_')[1]
    for record in SeqIO.parse(open(MSA, 'rU'), "fasta") :           
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
    test=TESTexist('/Users/steven/Desktop/Stage2016/SpecificFams/new/epicsMat/raxml/'+famname+'_mat_'+epicsTYPE)
    if test==False:
        return False
    f=open('/Users/steven/Desktop/Stage2016/SpecificFams/new/epicsMat/raxml/'+famname+'_mat_'+epicsTYPE)
    epicsmat=f.read().splitlines()
    for i in range(len(epicsmat)) :
        epicsmat[i]=epicsmat[i].split()
    f.close()
    return epicsmat
    
    
    
def get_epicsmatpval(famname,epicsTYPE):
    sumrow=[]
    test=TESTexist('/Users/steven/Desktop/Stage2016/SpecificFams/new/epicsMat/raxml/'+famname+'_matpval_'+epicsTYPE)
    if test==False:
        return False    
    f=open('/Users/steven/Desktop/Stage2016/SpecificFams/new/epicsMat/raxml/'+famname+'_matpval_'+epicsTYPE)
    epicsmat=f.read().splitlines()
    for i in range(len(epicsmat)) :
        epicsmat[i]=epicsmat[i].split()
        row = copy.deepcopy(epicsmat[i])
        row = map(float,row)
        sumrow.append(sum(row))
    f.close()
    return epicsmat,sumrow    



def parseSTRIDEfiles(STRIDEfile_list):
    STRIDE_dict={}
    stride_list=[]
    f=open(STRIDEfile_list,'r')
    stride_list=f.read().splitlines()
    f.close()
    for i in stride_list :
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
                AA=i[1]
    return structAA,enfouiAA,AA



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
    f=open('/Users/steven/Desktop/Stage2016/SpecificFams/new/RES/raxml/ '+familyname,'w')
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
    

     
#def do_distance_analysis(fam,seqs,threshold,typemat):
#    epicsmat,maxmat=get_epicsmat(fam,typemat)
#    epicspval,maxpval=get_epicsmatpval(fam,typemat)
#    family={}
#    for seq in seqs :
#        if 'cath' in seq.id :
#            goodAAs=[]
#            seqname=seq.id.split('|')[2].split('/')[0]
#            print seqname
#            print seq.seq
#            distmat=get_distmat(seqname)
#            print len(distmat)
#            for posAA1 in range(len(seq.seq)) :
#                for posAA2 in range(posAA1,len(seq.seq)):
#                    structAA1,structAA2=get_2nd_struct(seqname,posAA1,posAA2)
#                    if int(epicsmat[posAA1][posAA2]) > int(maxmat)/2 : #and distmat[posAA1][posAA2] < threshold  and epicspval[posAA1][posAA2] < threshold
#                       goodAAs.append(((seq.seq[posAA1],posAA1),(seq.seq[posAA2],posAA2),epicsmat[posAA1][posAA2],epicspval[posAA1][posAA2],(structAA1,structAA2)))
#                       print (seq.seq[posAA1],posAA1,structAA1,seq.seq[posAA2],posAA2,structAA2,epicsmat[posAA1][posAA2],epicspval[posAA1][posAA2],distmat[posAA1][posAA2])
#                       family[seqname]=goodAAs                  
#    return family                     
                        



def do_len_trimming():
    f=open('/Users/steven/Desktop/Stage2016/listsmall','r')
    listall=f.read().splitlines()
    f.close()
    listsmall=[]
    for i in range(len(listall)):
        listall[i]=listall[i].split()
        if int(listall[i][0]) <= 100 :
            listall[i][0]=int(listall[i][0])
            listsmall.append(listall[i])
    return listsmall
                


                            
def do_distance_analysis_bis(fam,seqs,atomseqs,cathdict,list_distmat):
    epicsmatCOOC=get_epicsmat(fam,'cooc')                 
    epicspvalCOOC,sumrowCOOC=get_epicsmatpval(fam,'cooc')
        
    #epicsmatCH=get_epicsmat(fam,'ch')                    
    #epicspvalCH,sumrowCH=get_epicsmatpval(fam,'ch')
    #if epicsmatCH==False or epicspvalCH==False :
        #return False
        
    #epicsmatCHCOOC=get_epicsmat(fam,'chcooc')                    
    #epicspvalCHCOOC,sumrowCHCOOC=get_epicsmatpval(fam,'chcooc')
    #if epicsmatCHCOOC==False or epicspvalCHCOOC==False :
        #return False
        
    goodAAsCOOC=[]
    #goodAAsCH=[]
    #goodAAsCHCOOC=[]
    
    Sequence={}
    for seq in seqs :
        if 'cath' in seq.id :
            goodAAsCOOC=[]
            #goodAAsCH=[]
            #goodAAsCHCOOC=[]
            print seq.id
            seqname=seq.id.split('_')[4]
            print seqname
            #Sequence[seqname]={'COOC':{},'CH':{},'CHCOOC':{}}
            Sequence[seqname]={'COOC':{}}
            alignedseq=AtC.do_mafft_addfull_keeplength(atomseqs[seqname],cathdict,list_distmat)
            print seqname
            print alignedseq.seq
            print len(seq.seq)
            distmat=get_distmat(seqname)
            if distmat==False :
                return False
            print 'get_distmat DONE'
            distmat=AtC.make_fullmat(distmat)
            print 'make_fullmat DONE'
            distmat=AtC.make_symetricalmat(distmat)
            print 'make_symetricalmat DONE'
            distmat=AtC.reformat_mat(distmat,alignedseq.seq)
            print 'reformat_mat DONE'
            AtC.write_Mat_file(distmat,seqname)
            print 'GETTING GOOD AAs'
            for posAA1 in range(len(seq.seq)):
                for posAA2 in range(posAA1,len(seq.seq)):
                    if sumrowCOOC[posAA1]==len(seq.seq) :
                        continue
                    gaps1=getnbGAPS(alignedseq.seq,posAA1)
                    gaps2=getnbGAPS(alignedseq.seq,posAA2)
                    if alignedseq.seq[posAA1]!='-':
                        structAA1,enfoui1,AA1stride=get_2nd_struct_bis(seqname,posAA1-gaps1+1)
                    else :
                        structAA1=0
                        enfoui1=0
                        AA1stride='-'
                    if alignedseq.seq[posAA2]!='-':
                        structAA2,enfoui2,AA2stride=get_2nd_struct_bis(seqname,posAA2-gaps2+1)
                    else :
                        structAA2=0
                        enfoui2=0
                        AA2stride='-'
                    if float(epicspvalCOOC[posAA1][posAA2]) < 0.01 and posAA1!=posAA2:
                        #print 'aligned seq pos1:'+ str(alignedseq.seq[posAA1])
                        #print 'stride pos1 :' + str(AA1stride)
                        #print 'aligned seq pos2:'+ str(alignedseq.seq[posAA2])
                        #print 'stride pos2 :' + str(AA2stride)
                        goodAAsCOOC.append((seq.seq[posAA1],posAA1,structAA1,enfoui1,seq.seq[posAA2],posAA2,structAA2,enfoui2,distmat[posAA1][posAA2],epicsmatCOOC[posAA1][posAA2],epicspvalCOOC[posAA1][posAA2]))
                    #if float(epicspvalCH[posAA1][posAA2]) < 0.05 and posAA1!=posAA2:
                        #goodAAsCH.append((seq.seq[posAA1],posAA1,structAA1,enfoui1,seq.seq[posAA2],posAA2,structAA2,enfoui2,distmat[posAA1][posAA2],epicsmatCH[posAA1][posAA2],epicspvalCH[posAA1][posAA2]))
                    #if float(epicspvalCHCOOC[posAA1][posAA2]) < 0.05 and posAA1!=posAA2:    
                        #goodAAsCHCOOC.append((seq.seq[posAA1],posAA1,structAA1,enfoui1,seq.seq[posAA2],posAA2,structAA2,enfoui2,distmat[posAA1][posAA2],epicsmatCHCOOC[posAA1][posAA2],epicspvalCHCOOC[posAA1][posAA2]))
            print 'GETTING GOOD AAs DONE'
            Sequence[seqname]['COOC']=goodAAsCOOC
            #Sequence[seqname]['CH']=goodAAsCH
            #Sequence[seqname]['CHCOOC']=goodAAsCHCOOC
            os.remove('/Users/steven/Desktop/Stage2016/dist_struct_mat_ATOM_new/'+seqname)
    return Sequence




def MAIN():
    #listfams=do_len_trimming()
    #listfamssorted=sorted(listfams, key=operator.itemgetter(0))
    #listfams=glob.glob('/Users/steven/Desktop/Stage2016/DATA/CATH_DATA/SequenceBySuperfamily/COMBS_CONSENSUS_TRIM/NR/MSA_MA/MSA_*')
    listfams=glob.glob('/Users/steven/Desktop/Stage2016/SpecificFams/new/MSAcleaned/MSA_*')[1:]
    list_distmat=os.listdir('/Users/steven/Desktop/Stage2016/dist_struct_Mat/')[1:]
    cathdict=AtC.parse_CATHlist('/Users/steven/Desktop/Stage2016/DATA/CATH_DATA/CathDomainList.S100.txt')
    atomseqs=getATOMS()
    families={}
    print listfams
    for i in listfams :
        #NAME=i.split('/')[7].split('_')[1]
        #PATH='/Users/steven/Desktop/Stage2016/Stage2016/SpecificFams/new/'
        fam,seqs=parseMSAfile(i)
        print fam
        families[fam]=do_distance_analysis_bis(fam,seqs,atomseqs,cathdict,list_distmat)
        if families[fam]==False:
            print 'END'
            continue
        writeRES(families[fam],fam)

MAIN()




       
#fam,seqs=parseMSAfile('/Users/steven/Desktop/Stage2016/DATA/CATH_DATA/SequenceBySuperfamily/COMBS_CONSENSUS_TRIM/NR/MSA_MA/','MSA_MA_S100.1.10.8.270_nC.lsf.nrpdb.rep.fasta.txt')
#stride=parseSTRIDEfiles('/Users/steven/Desktop/Stage2016/listStride')
#do_distance_analysis(fam,seqs,7,'COOC')
       
       
#get_positions(seqs[0].id.split('|')[2].split('/')[0])   
#f=open('/Users/steven/Desktop/Stage2016/listallCATH','r')
#listallCATH=f.read().splitlines()
#f.close() 
#do_distmatrix(listallCATH)