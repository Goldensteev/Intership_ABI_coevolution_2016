# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 15:20:22 2016

@author: steven
"""

import os


#listfiles=os.listdir('/Users/steven/Desktop/Stage2016/epicsTrees/')[1:]
#print len(listfiles)
#listall=[]
#
#l=os.listdir('/Users/steven/Desktop/Stage2016/DATA/CATH_DATA/SequenceBySuperfamily/COMBS_CONSENSUS_TRIM/NR/MSA_MA/')
#
#
#for i in l :
#    if 'renamedseqs_' in i :
#        listall.append(i.split('_')[3])     
#
#print len(listall)
#
#f=open('/Users/steven/Desktop/Stage2016/notin.txt','w')
#c=0
#for i in listall :
#    if i+'_epicsTree.txt' not in listfiles :
#
#        f.write('renamedseqs_MSA_MA'+i+'_nC.lsf.nrpdb.rep.fasta.txt\n')
#        c+=1
#f.close()

#print listall
#
#l=f.readlines()
#f.close()

l=os.listdir('/Users/steven/Desktop/Stage2016/DATA/CATH_DATA/SequenceBySuperfamily/COMBS_CONSENSUS_TRIM/epicsTrees/')[1:]
n=[]
for i in l :
    n.append(i[:-4])
    os.system('mv /Users/steven/Desktop/Stage2016/DATA/CATH_DATA/SequenceBySuperfamily/COMBS_CONSENSUS_TRIM/epicsTrees/'+i+' /Users/steven/Desktop/Stage2016/epicsTrees/'+i[:-4])
