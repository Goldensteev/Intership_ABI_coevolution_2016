# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 13:55:48 2016

@author: steven
"""

import os

listfiles=os.listdir('./DATA/')
for i in listfiles :
    if 'renamedseqs_' in i :
        print i
        if not os.path.exists('./FASTML_RUNS/'):
            os.makedir('./FASTML_RUNS/')
        os.mkdir('./FASTML_RUNS/'+i.split('_')[3]+'/')
        os.system('perl ./FastML.v3.1/www/fastml/FastML_Wrapper.pl --MSA_File ./DATA/'+i+' --seqType aa --outDir ./FASTML_RUNS/'+i.split('_')[3]+'/')
    