# -*- coding: utf-8 -*-
"""
Created on Tue May 24 10:41:48 2016

@author: steven
"""

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


change_name_seq('MSA_2.40.180.10.txt','/Users/steven/Desktop/Stage2016/SpecificFams/new/2.40.180.10/')