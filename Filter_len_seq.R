args = commandArgs(trailingOnly=TRUE)

data=read.table(args[1],sep = ' ')
data_m=mean(as.numeric(data$V2))
data_sd=sd(as.numeric(data$V2))
good_data=c()
if(length(data)!=1){
  for (i in (1:length(data$V1)))
  {
    if (data_m+data_sd > as.numeric(data$V2[i]) && as.numeric(data$V2[i]) > data_m-data_sd ) {
      good_data[length(good_data)+1]=as.character(data$V1[i])
    }
  }
  print('Y')
write(good_data, file = args[2])
}else{
  good_data=as.character(data$V1)
  write(good_data, file = args[2])
  print('N')
}