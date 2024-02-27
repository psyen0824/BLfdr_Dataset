# R Code for the Frequentist Bivariate Linear Mixed-Effects Model

# Import Data
ImagingData <- read.table ( "Github/Deidentified Data.txt", sep="", header = TRUE)

######################################################################## DTI Data Processing
ImagingData$Value[ImagingData$Measure=="DTI"] <- ImagingData$Value[ImagingData$Measure=="DTI"]^(1/3)

vectomat <- function(vec) matrix(vec[c(1,2,2,3)], nrow = 2)
mattovec <- function(mat) as.vector(mat)[c(1, 2, 4)]

numRegion <- length(unique(ImagingData$ROI_1)) + 1
numRegComb <- nrow(ImagingData)/numSubj/2
numSubj <- length(unique(ImagingData$ID))
numSubj0 <- as.numeric((table(ImagingData$Group)/3741/2)[1])
numSubj1 <- as.numeric((table(ImagingData$Group)/3741/2)[2])
GroupSubj <- c ( rep ( 0 , numSubj0 ) , rep ( 1 , numSubj1 ) )
R1R2 <- unique(ImagingData$ROI_1_ROI_2)

# Initializing Parameters
Beta <- matrix(0,numRegComb,4)
Beta.var <-matrix(0,numRegComb,16)
Beta.se <- matrix(0,numRegComb,4)
Random <- matrix(0,(numSubj0+numSubj1),2)
Sigma <- matrix(0,numRegComb,6)
Sigma.ran <- c(0.005,0,0.02)
sigma.ran.y <- array(0,c(2,2,(numSubj0+numSubj1)))
for(iSubj in 1:(numSubj0+numSubj1)) sigma.ran.y[,,iSubj] <- diag(c(0.005,0.02))

# Convergence Criteria
tol  <- 1e-5 
diff <- rep (FALSE, numRegComb*4+(numSubj0 + numSubj1)*2+numRegComb*6+4)
iter <- 1

Y <- as.numeric(ImagingData$Value)                      
X <- cbind(as.numeric(ImagingData$Measure=="fMRI" & ImagingData$Group==0), 
           as.numeric(ImagingData$Measure=="DTI"  & ImagingData$Group==0),
           as.numeric(ImagingData$Measure=="fMRI" & ImagingData$Group==1), 
           as.numeric(ImagingData$Measure=="DTI"  & ImagingData$Group==1)) 
Z <- cbind (as.numeric(factor(ImagingData$Measure)) - 1 ,2 - as.numeric (factor(ImagingData$Measure)))

while(iter < 100 && sum(diff) < length(diff)){
  print(paste("Iteration",iter))
  # Update ML Estimates of Beta,Sigma,and Sigma.ran
  sigma.ran.mat <-t(Random)%*%Random/(numSubj0+numSubj1)+apply(sigma.ran.y,c(1,2),mean)
  Sigma.ran.new <- mattovec(sigma.ran.mat)
  Beta.new  <- matrix(0,numRegComb,4)
  Sigma.new <- matrix(0,numRegComb,6)
  sigma0.inv <- matrix(0,2*numRegComb,2*numRegComb)
  sigma1.inv <- matrix(0,2*numRegComb,2*numRegComb)
  
  for(iRegComb in 1:numRegComb){
    temp <- ImagingData[ImagingData$ROI_1_ROI_2==R1R2[iRegComb],]
    y <- Y[ImagingData$ROI_1_ROI_2==R1R2[iRegComb]]
    x <- X[ImagingData$ROI_1_ROI_2==R1R2[iRegComb],]
    z <- Z[ImagingData$ROI_1_ROI_2==R1R2[iRegComb],]
    error0 <- matrix((y-x%*%Beta[iRegComb,]-as.numeric(Random))[temp$Group==0],ncol=2)
    error1 <- matrix((y-x%*%Beta[iRegComb,]-as.numeric(Random))[temp$Group==1],ncol=2)
    
    sigma.mat0 <- (t(error0)%*%error0+apply(sigma.ran.y[,,GroupSubj==0],c(1,2),sum))/numSubj0
    sigma.mat1 <- (t(error1)%*%error1+apply(sigma.ran.y[,,GroupSubj==1],c(1,2),sum))/numSubj1
    
    Sigma.new[iRegComb,]<-c(mattovec(sigma.mat0),mattovec(sigma.mat1))
    
    sigma0 <- c(Sigma.new[iRegComb,3],-Sigma.new[iRegComb,2],Sigma.new[iRegComb,1]) /(Sigma.new[iRegComb,1]*Sigma.new[iRegComb,3]-Sigma.new[iRegComb,2]^2)
    sigma1 <- c(Sigma.new[iRegComb,6],-Sigma.new[iRegComb,5],Sigma.new[iRegComb,4]) /(Sigma.new[iRegComb,4]*Sigma.new[iRegComb,6]-Sigma.new[iRegComb,5]^2)
   
    # Define Error Covariance Matrix
    sigma0.inv[iRegComb,iRegComb] <- sigma0[1]
    sigma0.inv[iRegComb+numRegComb,iRegComb+numRegComb]<- sigma0[3]
    sigma0.inv[iRegComb,iRegComb+numRegComb]<- sigma0[2]
    sigma0.inv[iRegComb+numRegComb,iRegComb]<- sigma0[2]
    sigma1.inv[iRegComb,iRegComb]<- sigma1[1]
    sigma1.inv[iRegComb+numRegComb,iRegComb+numRegComb]<- sigma1[3]
    sigma1.inv[iRegComb,iRegComb+numRegComb]<- sigma1[2]
    sigma1.inv[iRegComb+numRegComb,iRegComb]<- sigma1[2]
    sigma.e <- matrix(0,(numSubj0+numSubj1)*2,(numSubj0+numSubj1)*2)
    
    for(iSubj in 1:(numSubj0+numSubj1)){
      sigma.e[iSubj,iSubj]<-Sigma.new[iRegComb,1]*(1-GroupSubj[iSubj]) + Sigma.new[iRegComb,4]*GroupSubj[iSubj]
      sigma.e[((numSubj0+numSubj1)+iSubj),((numSubj0+numSubj1)+iSubj)]<- Sigma.new[iRegComb,3]*(1-GroupSubj[iSubj]) + Sigma.new[iRegComb,6]*GroupSubj[iSubj]
      sigma.e[iSubj,((numSubj0+numSubj1)+iSubj)]<- Sigma.new[iRegComb,2]*(1-GroupSubj[iSubj]) +Sigma.new[iRegComb,5]*GroupSubj[iSubj]
      sigma.e[((numSubj0+numSubj1)+iSubj),iSubj]<- sigma.e[iSubj,((numSubj0+numSubj1)+iSubj)]
    }
    
    Beta.new[iRegComb,]<- solve(t(x)%*%sigma.e%*%x) %*% t(x) %*% sigma.e %*%(y-as.numeric(Random))
    
    Beta.var[iRegComb,] <- as.numeric(solve(t(x)%*%solve(sigma.e + z%*%sigma.ran.mat %*%t(z)) %*%x))
    
    Beta.se[iRegComb,]<- sqrt(diag(solve(t(x)%*%solve(sigma.e + z%*%sigma.ran.mat %*% t(z)) %*% x)))
  }
  
  #Update EB Estimates and Variance of Random Effects
  Random.new <- matrix(0,(numSubj0+numSubj1),2)
  
  for(iSubj in 1:(numSubj0+numSubj1)){
    temp <- ImagingData[ImagingData$ID==iSubj,]
    y <- Y[ImagingData$ID==iSubj]
    x <- X[ImagingData$ID==iSubj,]
    z <- Z[ImagingData$ID==iSubj,]
    
    if(unique(temp$Group)==0) r <- vectomat(Sigma.ran.new)%*%solve(vectomat(Sigma.ran.new)+solve(t(z)%*%sigma0.inv%*%z)) else r <- vectomat(Sigma.ran.new)%*%solve(vectomat(Sigma.ran.new)+solve(t(z)%*%sigma1.inv%*%z))
    
    Random.new[iSubj,]<- r%*%solve(t(z)%*%z)%*%t(z)%*%(-diag(x%*%t(Beta.new[rep(1:numRegComb,2),])))
    sigma.ran.y[,,iSubj]<-(diag(2)-r)%*%vectomat(Sigma.ran.new)
    
  }
  
  diff<-abs(c(as.numeric(Beta.new-Beta),as.numeric(Random.new-Random),
              as.numeric(Sigma.new-Sigma),as.numeric(Sigma.ran.new-Sigma.ran)))<tol
  print(paste(sum(diff),"out of",length(diff),"parameters converged"))
  Beta <- Beta.new
  Random <- Random.new
  Sigma <-Sigma.new
  Sigma.ran <-Sigma.ran.new
  iter <- iter + 1
}
