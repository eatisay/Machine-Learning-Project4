dataSet<-read.csv("hw04_data_set.csv")
#Splitting the data
trainSet<-dataSet[c(1:150),]
testSet<-dataSet[c(151:272),]
dataSetW<-dataSet$waiting
dataSetE<-dataSet$eruptions
trainW<-trainSet$waiting
trainE<-trainSet$eruptions
testW<-testSet$waiting
testE<-testSet$eruptions
#Parameters to use later on
binWidth <- 0.37
originParam<-1.5
minimumValue<-min(dataSetE)
maximumValue<-max(dataSetE)
#Data interval to use in the running mean smoother and kernel smoother
dataInterval <- seq(from = originParam, to = maximumValue, by = 0.01)
#REGRESSOGRAM
#Data interval from the min and max values to use in regressogram
leftBorders <- seq(from = originParam, to = maximumValue, by = binWidth)
rightBorders <- seq(from = originParam + binWidth, to = maximumValue+ binWidth, by = binWidth)
#Prediction is calcuated as means of the values in the same bin. 
gHeadR <- sapply(1:length(leftBorders), function(x) {
  return(mean(trainW[leftBorders[x] < trainE & trainE <= rightBorders[x]]))
}
)
#Plotting
plot(x=trainE,y=trainW,xlab="Eruption time (min)", ylab="Waiting time to next eruption (min)",
     type = "p", main = sprintf("h = %g", binWidth),col = "purple",las = 1,pch = 19
   )
legend("topleft",legend=c("training", "test"),fill = c("purple","red"), cex=0.8,pt.cex = 0.001)
points(x=testE,y=testW, col = "red",las = 1,pch = 19)
for (b in 1:length(leftBorders)) {
  lines(c(leftBorders[b], rightBorders[b]), c(gHeadR[b], gHeadR[b]), lwd = 2, col = "black")
  if (b < length(leftBorders)) {
    lines(c(rightBorders[b], rightBorders[b]), c(gHeadR[b], gHeadR[b + 1]), lwd = 2, col = "black") 
  }
}
#Calculating the RMSE
testLength<-length(testE)
i<-1
nominator<-c()
while (i<=length(testE)) {
  #Getting the corresponding bin
  bin<-(testE[i]-minimumValue) / binWidth
  #Calculating the deviance
  nom<-(testW[i]-gHeadR[bin])^2
  nominator<-cbind(nominator,nom)
  i<-i+1
}
RMSEforR<-sqrt(sum(nominator)/testLength)
print(paste("Regressogram => RMSE is", RMSEforR, "when h is", binWidth))
#RUNNING MEAN SMOOTHER
#w function for having a symmteric bin around a point in data inerval
funcW<-function(x){if(abs(x)<1)return(1) else return(0)}
funcW<-function(x){if(abs(x)<1/2)return(1) else return(0)}
#Sum of the values from that bin
funcRMSNominator<-function(x,l){
  res<-c()
  i=1
  while (i<=l) {
    v <- (x - trainE[i]) / binWidth
    weight <- funcW(v)
    res<-cbind(res,(weight*trainW[i]))
    i<-i+1
  }
  return(sum(res))
}
#Number of values in that bin
funcRMSDenominator<-function(x,l){
  res<-c()
  i=1
  while (i<=l) {
    v <- (x - trainE[i]) / binWidth
    weight <- funcW(v)
    res<-cbind(res,weight)
    i<-i+1
  }
  return(sum(res))
}
#Predicted values for bins from the data interval
gHeadRMS <- sapply(dataInterval, function(x) {
  nominator <- funcRMSNominator(x,length(trainE))
  denominator <- funcRMSDenominator(x,length(trainE))
  return (nominator / denominator)
})
#Plotting
plot(x=trainE,y=trainW,xlab="Eruption time (min)", ylab="Waiting time to next eruption (min)",
     type = "p", main = sprintf("h = %g", binWidth),col = "purple",las = 1,pch = 19
)
legend("topleft",legend=c("training", "test"),fill = c("purple","red"), cex=0.8,pt.cex = 0.001)
points(x=testE,y=testW, col = "red",las = 1,pch = 19)
for (b in 1:length(dataInterval)) {
  lines(c(dataInterval[b], dataInterval[b+1]), c(gHeadRMS[b], gHeadRMS[b]), lwd = 2, col = "black")
  if (b < length(dataInterval)) {
    lines(c(dataInterval[b], dataInterval[b+1]), c(gHeadRMS[b], gHeadRMS[b + 1]), lwd = 2, col = "black") 
  }
}
#Calculating the RMSE
i<-1
nominator<-c()
while (i<=length(testE)) {
  #Calculating the corresponding bin
  inter<-(testE[i]-minimumValue) / 0.01 #Denominator varies as incrementation in data interval's sequence
  nom<-(testW[i]-gHeadRMS[inter])^2
  nominator<-cbind(nominator,nom)
  i<-i+1
}
RMSEforRMS<-sqrt(sum(nominator)/testLength)
print(paste("Running Mean Smoother => RMSE is", RMSEforRMS, "when h is", binWidth))
#KERNEL SMOOTHER
#Calculating the Gaussian Kernel from the lecture notes
funcKernel<- function(v){(1/sqrt(2*pi))*exp(-v^2/2)}
#Nominator of the prediction function., resulting the knn smoothing
funcKernelNominator<-function(x,l){
  res<-c()
  i=1
    while (i<=l) {
      v <- (x - trainE[i]) / binWidth
      kernel <- funcKernel(v)
      res<-cbind(res,(kernel*trainW[i]))
      i<-i+1
    }
  return(sum(res))
}
#Denominator for the prediction.
funcKernelDenominator<-function(x,l){
  res<-c()
  i=1
  while (i<=l) {
    v <- (x - trainE[i]) / binWidth
    kernel <- funcKernel(v)
    res<-cbind(res,kernel)
    i<-i+1
  }
  return(sum(res))
}
#Prediction for the intervals using Kernel smoothining
gHeadKS <- sapply(dataInterval, function(x) {
  nominator <- funcKernelNominator(x,length(trainE))
  denominator <- funcKernelDenominator(x,length(trainE))
  return (nominator / denominator)
})
#Plotting
plot(x=trainE,y=trainW,xlab="Eruption time (min)", ylab="Waiting time to next eruption (min)",
     type = "p", main = sprintf("h = %g", binWidth),col = "purple",las = 1,pch = 19
)
legend("topleft",legend=c("training", "test"),fill = c("purple","red"), cex=0.8,pt.cex = 0.001)
points(x=testE,y=testW, col = "red",las = 1,pch = 19)
for (b in 1:length(dataInterval)) {
  lines(c(dataInterval[b], dataInterval[b+1]), c(gHeadKS[b], gHeadKS[b]), lwd = 2, col = "black")
  if (b < length(dataInterval)) {
    lines(c(dataInterval[b], dataInterval[b+1]), c(gHeadKS[b], gHeadKS[b + 1]), lwd = 2, col = "black") 
  }
}
#Calculating the RMSE
i<-1
nominator<-c()
while (i<=length(testE)) {
  #Getting the corresponding interval.
  inter<-(testE[i]-minimumValue) / 0.01 #Again, this denominator changes as the data interval
  nom<-(testW[i]-gHeadKS[inter])^2
  nominator<-cbind(nominator,nom)
  i<-i+1
}
RMSEforKS<-sqrt(sum(nominator)/testLength)
print(paste("Kernel Smoother => RMSE is", RMSEforKS, "when h is", binWidth))

