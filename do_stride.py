# -*- coding: utf-8 -*-
"""
Created on Mon May  9 10:25:37 2016

@author: steven
"""

import os 

def do_stride(dompdb_list):
    for i in dompdb_list :
        print i
        #print '/Users/steven/Desktop/Stage2016/stride/./stride /Users/steven/Desktop/Stage2016/DATA/CATH_DATA/dompdb/'+i+'> /Users/steven/Desktop/Stage2016/DATA/CATH_DATA/Stride_dompdb_out/'+i
        os.system('/Users/steven/Desktop/Stage2016/stride/./stride /Users/steven/Desktop/Stage2016/DATA/CATH_DATA/dompdb/'+i+'>/Users/steven/Desktop/Stage2016/DATA/CATH_DATA/Stride_dompdb_out/'+i)

        
def read_DOMPDBlist(DOMPDBfile):
    dompdb_list=[]     
    f=open(DOMPDBfile,'r')
    dompdb_list=f.read().splitlines()
    f.close()
    return dompdb_list
    
    
def MAIN():
    dompdb_list=read_DOMPDBlist('/Users/steven/Desktop/Stage2016/listdompdb')
    do_stride(dompdb_list)

MAIN()     
