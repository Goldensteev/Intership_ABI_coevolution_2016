library(lattice)

data=read.table(file='/Users/steven/Desktop/Stage2016/epicsMAT/CH+COOC/S100.1.10.8.20_matPVAL')
d=as.matrix(data)
levelplot(d)
 
