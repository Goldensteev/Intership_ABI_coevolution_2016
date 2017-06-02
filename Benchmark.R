data1=read.table(file = '/Users/steven/Desktop/Stage2016/FASTMLbenchmarking/BenchmarkFASTMLsource.Lseq=350.txt',header = TRUE)
data2=read.table(file = '/Users/steven/Desktop/Stage2016/FASTMLbenchmarking/BenchmarkFASTMLmytimeLseq=721.txt',header = TRUE)
log.data1=data1
log.data2=data2
log.data1$NBmin=log(log.data1$NBmin)
model1 = lm(log.data1$NBmin ~ log.data1$NBseq, data1 = log.data1)
plot(data1)
lines(data1$NBseq,exp(fitted(model1)))


log.data2$NBmin=log(log.data2$NBmin)
model2 = lm(log.data2$NBmin ~ log.data2$NBseq, data2 = log.data2)
plot(data2)
lines(data2$NBseq,exp(fitted(model2)))
