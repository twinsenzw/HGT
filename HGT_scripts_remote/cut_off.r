args = commandArgs(trailingOnly=TRUE)
distance <- scan(args[1], what="", sep="\n")
distance.numeric <- as.numeric(distance)

max = 50

distance.nout <- c()
for(i in 1:length(distance.numeric))
{
  if(distance.numeric[i]<=max)
  {
    distance.nout <- append(distance.nout,distance.numeric[i])
  }
}

cutoff <- qnorm(0.05,mean=mean(distance.nout),sd=sd(distance.nout))
write(cutoff,args[2])


