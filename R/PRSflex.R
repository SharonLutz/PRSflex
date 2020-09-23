PRSflex <-
function(nSim=50, nTrain=2500, nTest = 2500, varY=1, sigl=5e-8, betaS = seq(0,0.1,by=0.05),
                   nAdd=25, nDom=25, nCod=25, nRec=25,
                   MAFadd = 0.05, MAFdom = 0.05, MAFcod = 0.05, MAFrec = 0.05,
                   plot.pdf=TRUE, SEED=1){

  library(caret) # for preProcess
  #library(stats) #for AIC
  library(HDeconometrics) # for lasso
  
  set.seed(SEED) # set the seed
  nSNP<-nAdd+nDom+nCod+nRec #total number of snps
  n<- nTrain+nTest #total sample size    
  minMAF <- min(MAFadd, MAFdom, MAFcod, MAFrec)
  
  ##################################################################
  # Error check
  ##################################################################
  if(varY <= 0){stop("Error: varY must be greater than 0")}
  if(MAFadd >= 1 | MAFadd <= 0){stop("Error: MAFadd must be greater than 0 and less than 1")}
  if(MAFdom >= 1 | MAFdom <= 0){stop("Error: MAFdom must be greater than 0 and less than 1")}
  if(MAFcod >= 1 | MAFcod <= 0){stop("Error: MAFcod must be greater than 0 and less than 1")}
  if(MAFrec >= 1 | MAFrec <= 0){stop("Error: MAFrec must be greater than 0 and less than 1")}
  if(nTrain+nTest < 40){stop("Error: not enough people, nTrain+nTest < 40")}
  if(nAdd+nDom+nCod+nRec < 40){stop("Error: not enough people, nAdd+nDom+nCod+nRec < 40")}
  if(minMAF*nTrain < 5){stop("Error: minimum MAF x nTrain must be greater than 10")}
  if(minMAF*nTest < 5){stop("Error: minimum MAF x nTest must greater than 10")}
  
  ##################################################################
  #store results
  ##################################################################
  matR<-matrix(0,nrow=length(betaS),ncol=18)
  colnames(matR)<-c("cLRmse","cLASSOmse","pLRmse","pLASSOmse","p2LRmse","p2LASSOmse",
                    "cLRadjR2","cLASSOadjR2","pLRadjR2","pLASSOadjR2","p2LRadjR2","p2LASSOadjR2",
                    "cLRp","cLASSOp","pLRp","pLASSOp","p2LRp","p2LASSOp")
  
  # c stands for classic and p stands for proposed
  # LR stands for linear regression 
  # Lasso is the other approach
  #then is gives the mse and adjusted R squared for all 4 methods
  
  ##################################################################
  # loop betaS
  ##################################################################
  for(bb in 1:length(betaS)){
    
    ##################################################################
    # loop sims
    ##################################################################
    for(ss in 1:nSim){
      #print every 10 simulations
      printS<-10
      if(floor(ss/printS)==ceiling(ss/printS)){print(paste(bb,"in",length(betaS),"and",ss,"in",nSim))}
      
      ##################################################################
      # Generate the data
      ##################################################################
      X<-c() # snps as 0,1,2 which is known
      Xtrue<-c() #true genetic coding which is unknown
      Xind<-c() #indicator variables for X
      
      #additive ggenetic coding
      if(nAdd>0){
        for(aa in 1:nAdd){
          snp<-rbinom(n,2,MAFadd)
          X<-cbind(X,snp)
          Xtrue<-cbind(Xtrue,snp)
          
          #indicator variables only
          snp1<-rep(0,n)
          snp2<-rep(0,n)
          snp1[snp==1]<-1
          snp2[snp==2]<-1
          snpv<-cbind(snp1,snp2)
          Xind<-cbind(Xind,snpv)
        }}
      #dominant 1&2 vs 0
      if(nDom>0){
        for(dd in 1:nDom){
          snp<-rbinom(n,2,MAFdom)
          X<-cbind(X,snp)
          snp2<-snp
          snp2[snp==2]<-1
          Xtrue<-cbind(Xtrue,snp2)
          
          #indicator variables only
          snp1<-rep(0,n)
          snp2<-rep(0,n)
          snp1[snp==1]<-1
          snp2[snp==2]<-1
          snpv<-cbind(snp1,snp2)
          Xind<-cbind(Xind,snpv)
        }}
      #codominant
      if(nCod>0){
        for(cc in 1:nCod){
          snp<-rbinom(n,2,MAFcod)
          X<-cbind(X,snp)
          snp1<-rep(0,n)
          snp2<-rep(0,n)
          snp1[snp==1]<-1
          snp2[snp==2]<-1
          snpv<-cbind(snp1,snp2)
          Xtrue<-cbind(Xtrue,snpv)
          
          #indicator variables only
          snp1<-rep(0,n)
          snp2<-rep(0,n)
          snp1[snp==1]<-1
          snp2[snp==2]<-1
          snpv<-cbind(snp1,snp2)
          Xind<-cbind(Xind,snpv)
        }}
      #recessive 0&1 vs 2
      if(nRec>0){
        for(rr in 1:nRec){
          snp<-rbinom(n,2,MAFrec)
          X<-cbind(X,snp)
          snp2<-snp
          snp2[snp==1]<-0
          Xtrue<-cbind(Xtrue,snp2)
          
          #indicator variables only
          snp1<-rep(0,n)
          snp2<-rep(0,n)
          snp1[snp==1]<-1
          snp2[snp==2]<-1
          snpv<-cbind(snp1,snp2)
          Xind<-cbind(Xind,snpv)
        }}
      
      #colnames for X to use in LASSO
      colnames(X)<-c(paste0("x", 1:ncol(X)))
      xnamX<-colnames(X)
      
      #colnames for Xind to use in LASSO
      colnames(Xind)<-c(paste0(c("x1ind","x2ind"), rep(1:ncol(X),each=2)))
      xnamXind<-colnames(Xind)
      
      #generate y based on the true genetic coding which is unknown
      y<-rnorm(n,rowSums(Xtrue*betaS[bb]),sd=sqrt(varY))
      
      
      # generate flexible X based on AIC for indicators or continous
      Xflex<-c()
      for(xx in 1:ncol(X)){
        AICcont<-AIC(lm(y~X[,xx])) #continous model
        AICind<-AIC(lm(y~as.factor(X[,xx])))
        
        if(AICind<AICcont){
          # indicator variables for 0 and 2 with reference group of 1
          # indicator for 2 is dropped then indicator 0 becomes 0 vs 1,2 dominant model
          # indicator for 0 is dropped then indicator 2 becomes 2 vs 1,0 recessive model
          x0<-x2<-rep(0,length(X[,xx]))
          x0[X[,xx]==0]<-1
          x2[X[,xx]==2]<-1
          
          Xflex<-cbind(Xflex,x0,x2)
          colnames(Xflex)[(ncol(Xflex)-1):ncol(Xflex)]<-c(paste(c("i0x","i2x"),xx,sep=""))
        }else{
          # continous model, not indicator variables
          Xflex<-cbind(Xflex,X[,xx])
          colnames(Xflex)[ncol(Xflex)]<-c(paste("x",xx,sep=""))
        }
      }
      
      #combine X and y for data formatting 
      Xy<-cbind(X,y) 
      
      ##################################################################
      # classic method: create training and test sets
      ##################################################################
      indx<-c(1:nTrain)
      trainData = Xy[indx,]
      testData = Xy[-indx,] # Create the testData
      trainData<-data.frame(trainData)
      testData<-data.frame(testData)
      
      ##################################################################
      # classic method: Linear regression
      ##################################################################
      # Step 1 - fit linear regression on training data
      lr = lm(as.formula(paste("y ~ ", paste(xnamX, collapse= "+"))),data=trainData)
      
      # Step 2 - LR predicting and evaluating the model on testData data
      predLR<-predict(object=lr,newdata=testData)
      modelPLR<-summary(lm(y~predLR,data=testData))
      
      # store adjusted R squared and MSE  and p-value for association with PRS and y
      matR[bb,"cLRadjR2"]<-matR[bb,"cLRadjR2"]+modelPLR$adj.r
      if(nrow(modelPLR$coef)>1 & modelPLR$coef[nrow(modelPLR$coef),4]<sigl){matR[bb,"cLRp"]<-matR[bb,"cLRp"]+1} 
      matR[bb,"cLRmse"]<-matR[bb,"cLRmse"]+mean((predLR-testData$y)^2)
      
      ##################################################################
      # classic method: Lasso
      ##################################################################
      #scale the data
      pre_proc_val <- preProcess(trainData[,1:(ncol(Xy)-1)], method = c("center", "scale"))
      
      #use scaled data
      trainData[,1:(ncol(Xy)-1)] = predict(pre_proc_val, trainData[,1:(ncol(Xy)-1)])
      testData[,1:(ncol(Xy)-1)] = predict(pre_proc_val, testData[,1:(ncol(Xy)-1)])
      
      # need matrix form for glmnet
      trainData<-as.matrix(trainData)
      testData<-as.matrix(testData)
      
      # lasso model then predict - specify AIC or BIC
      lasso_model <- HDeconometrics::ic.glmnet(x=trainData[,1:(ncol(Xy)-1)], y=trainData[,"y"], crit="aic",alpha=1, standardize=TRUE)
      predictions_test <- predict(lasso_model, newdata = testData[,1:(ncol(Xy)-1)])
      
      testData<-data.frame(testData)
      modelPLA<-summary(lm(y~predictions_test,data=testData))
      
      matR[bb,"cLASSOadjR2"]<-matR[bb,"cLASSOadjR2"]+modelPLA$adj.r
      if(nrow(modelPLA$coef)>1 & modelPLA$coef[nrow(modelPLA$coef),4]<sigl){matR[bb,"cLASSOp"]<-matR[bb,"cLASSOp"]+1}
      matR[bb,"cLASSOmse"]<-matR[bb,"cLASSOmse"]+mean((predictions_test-testData$y)^2)
      
      ##################################################################
      # proposed method: train and test sets
      ##################################################################
      #make sure there are not too many 0 or 1 for the indicators
      # ie get rid if indicator variables with very little variablity
      
      indx<-c(1:nTrain)
      Xflex2<-c()
      xnamXflex<-c()
      
      
      for(tt in 1:ncol(Xflex)){
        if(var(Xflex[indx, tt])>0.001 & var(Xflex[-indx, tt])>0.001){
          Xflex2<-cbind(Xflex2,Xflex[,tt])
          xnamXflex<-c(xnamXflex,colnames(Xflex)[tt])
        }
      }
      colnames(Xflex2) <- xnamXflex
      Xflexy<-cbind(Xflex2,y) 
      
      
      trainData2 = Xflexy[indx,]
      testData2 = Xflexy[-indx,] # Create the testData2
      trainData2<-data.frame(trainData2)
      testData2<-data.frame(testData2)
      
      
      ##################################################################
      # proposed method: Linear Regression
      ##################################################################
      # Step 1 - fit linear regression on training Data2
      lr = lm(as.formula(paste("y ~ ", paste(xnamXflex, collapse= "+"))),data=trainData2)
      
      # Step 2 - LR predicting and evaluating the model on testData2 Data2
      predLR<-predict(object=lr,newdata=testData2)
      modelPLR<-summary(lm(y~predLR,data=testData2))
      
      # store adjusted R squared and MSE  and p-value for association with PRS and y
      matR[bb,"pLRadjR2"]<-matR[bb,"pLRadjR2"]+modelPLR$adj.r
      if(nrow(modelPLR$coef)>1 & modelPLR$coef[nrow(modelPLR$coef),4]<sigl){matR[bb,"pLRp"]<-matR[bb,"pLRp"]+1} 
      matR[bb,"pLRmse"]<-matR[bb,"pLRmse"]+mean((predLR-testData2$y)^2)
      
      ##################################################################
      # proposed method: Lasso
      ##################################################################
      #scale the data
      pre_proc_val <- preProcess(trainData2[,1:(ncol(trainData2)-1)], method = c("center", "scale"))
      
      # if(!is.null(trainData2$x27)){
      #  print(1)
      #   print(var(trainData2[,"x27"]))}
      
      #use scaled data
      trainData2[,1:(ncol(trainData2)-1)] = predict(pre_proc_val, trainData2[,1:(ncol(trainData2)-1)])
      testData2[,1:(ncol(testData2)-1)] = predict(pre_proc_val, testData2[,1:(ncol(testData2)-1)])
      
      # need matrix form for glmnet
      trainData2<-as.matrix(trainData2)
      testData2<-as.matrix(testData2)
      
      # lasso model then predict - specify AIC or BIC
      lasso_model <- HDeconometrics::ic.glmnet(x=trainData2[,1:(ncol(trainData2)-1)], y=trainData2[,"y"], crit="aic",alpha=1, standardize=TRUE)
      predictions_test <- predict(lasso_model, newdata = testData2[,1:(ncol(trainData2)-1)])
      
      testData2<-data.frame(testData2)
      modelPLA<-summary(lm(y~predictions_test,data=testData2))
      
      
      matR[bb,"pLASSOadjR2"]<-matR[bb,"pLASSOadjR2"]+modelPLA$adj.r
      if(nrow(modelPLA$coef)>1 & modelPLA$coef[nrow(modelPLA$coef),4]<sigl){matR[bb,"pLASSOp"]<-matR[bb,"pLASSOp"]+1}
      matR[bb,"pLASSOmse"]<-matR[bb,"pLASSOmse"]+mean((predictions_test-testData2$y)^2)
      
      ##################################################################
      # classic method: create training and test sets
      ##################################################################
      #make sure there are not too many 0 or 1 for the indicators
      # ie get rid if indicator variables with very little variablity
      
      # increase mean cutoff to .01 or do check on test and train make sure in both
      # get rid of mean stuff and look at variance - variance below 0.001 exclude (do on both train and test)
      indx<-c(1:nTrain)
      Xind2<-c()
      xnamXind<-c()
      
      for(tt in 1:ncol(Xind)){
        if(var(Xind[indx,tt])>0.001 & var(Xind[-indx,tt])>0.001){
          Xind2<-cbind(Xind2,Xind[,tt])
          xnamXind<-c(xnamXind,colnames(Xind)[tt])
        }
      }
      Xindy<-cbind(Xind,y) 
      
      trainData2 = Xindy[indx,c(xnamXind, "y")]
      testData2 = Xindy[-indx, c(xnamXind, "y")] # Create the testData2
      trainData2<-data.frame(trainData2)
      testData2<-data.frame(testData2)
      
      ##################################################################
      # proposed method: Linear Regression
      ##################################################################
      # Step 1 - fit linear regression on training Data2
      lr = lm(as.formula(paste("y ~ ", paste(xnamXind, collapse= "+"))),data=trainData2)
      
      # Step 2 - LR predicting and evaluating the model on testData2 Data2
      predLR<-predict(object=lr,newdata=testData2)
      modelp2LR<-summary(lm(y~predLR,data=testData2))
      
      # store adjusted R squared and MSE  and p-value for association with PRS and y
      matR[bb,"p2LRadjR2"]<-matR[bb,"p2LRadjR2"]+modelp2LR$adj.r
      if(nrow(modelp2LR$coef)>1 & modelp2LR$coef[nrow(modelp2LR$coef),4]<sigl){matR[bb,"p2LRp"]<-matR[bb,"p2LRp"]+1} 
      matR[bb,"p2LRmse"]<-matR[bb,"p2LRmse"]+mean((predLR-testData2$y)^2)
      
      ##################################################################
      # proposed method: Lasso
      ##################################################################
      #scale the data
      pre_proc_val <- preProcess(trainData2[,1:(ncol(trainData2)-1)], method = c("center", "scale"))
      
      
      #use scaled data
      trainData2[,1:(ncol(trainData2)-1)] = predict(pre_proc_val, trainData2[,1:(ncol(trainData2)-1)])
      testData2[,1:(ncol(testData2)-1)] = predict(pre_proc_val, testData2[,1:(ncol(testData2)-1)])
      
      # need matrix form for glmnet
      trainData2<-as.matrix(trainData2)
      testData2<-as.matrix(testData2)
      
      # lasso model then predict - specify AIC or BIC
      lasso_model <- HDeconometrics::ic.glmnet(x=trainData2[,1:(ncol(trainData2)-1)], y=trainData2[,"y"], crit="aic",alpha=1, standardize=TRUE)
      predictions_test <- predict(lasso_model, newdata = testData2[,1:(ncol(trainData2)-1)])
      
      testData2<-data.frame(testData2)
      modelPLA<-summary(lm(y~predictions_test,data=testData2))
      
      matR[bb,"p2LASSOadjR2"]<-matR[bb,"p2LASSOadjR2"]+modelPLA$adj.r
      if(nrow(modelPLA$coef)>1 & modelPLA$coef[nrow(modelPLA$coef),4]<sigl){matR[bb,"p2LASSOp"]<-matR[bb,"p2LASSOp"]+1}
      matR[bb,"p2LASSOmse"]<-matR[bb,"p2LASSOmse"]+mean((predictions_test-testData2$y)^2)
      
      
      ##################################################################
      #end ss loop simmulations
      ##################################################################
    }
    
    ##################################################################
    #end bb loop betaS
    ##################################################################
  }
  ##################################################################
  ##################################################################
  # Results
  ##################################################################
  
  #results- make sure to devide by the number of simulations
  res<-matR/nSim

    write.table(res, file =paste0("res_", "nA",nAdd,"nD",nDom,"nC",nCod,"nR",nRec,
                                "MAFa",MAFadd*100,"MAFd",MAFdom*100,"MAFc",MAFcod*100,"MAFr",MAFrec*100, 
                                ".txt"), row.names = F, col.names = T, quote=F)

  ##################################################################
  # plot adjusted R^2
  ##################################################################
  
  if(plot.pdf==TRUE){
  
  adjR2Vals<-c(res[,"cLRadjR2"],res[,"cLASSOadjR2"],res[,"pLRadjR2"],res[,"pLASSOadjR2"],res[,"p2LRadjR2"],res[,"p2LASSOadjR2"])
  
  pdf(paste("plotAdjR2","nA",nAdd,"nD",nDom,"nC",nCod,"nR",nRec,
            "MAFa",MAFadd*100,"MAFd",MAFdom*100,"MAFc",MAFcod*100,"MAFr",MAFrec*100,".pdf",sep=""))
  plot(-1,-1,xlim=c(min(betaS),max(betaS)),ylim=c(0,max(adjR2Vals)+sd(adjR2Vals)),ylab="Average Adjusted R^2",
       xlab="Genetic Effect Size",cex=1.9)
  lines(betaS,res[,"cLRadjR2"],col=1,lty=3,type="b", pch=2,lwd=2)
  lines(betaS,res[,"cLASSOadjR2"],col=3,lty=3,type="b", pch=3,lwd=2)
  lines(betaS,res[,"pLRadjR2"],col=4,lty=3,type="b", pch=4,lwd=2)
  lines(betaS,res[,"pLASSOadjR2"],col=6,lty=3,type="b", pch=5,lwd=2)
  lines(betaS,res[,"p2LRadjR2"],col=7,lty=3,type="b", pch=6,lwd=2)
  lines(betaS,res[,"p2LASSOadjR2"],col=8,lty=3,type="b", pch=7,lwd=2)
  legend("topleft", inset = 0.01,c("Classic PRS: regression","Classic PRS: LASSO",
                                   "Proposed PRS 1: regression","Proposed PRS 1: LASSO",
                                   "Proposed PRS 2: regression","Proposed PRS 2: LASSO"),col=c(1,3,4,6,7,8),pch=c(2:7),lty=3,lwd=2,cex=1)
  dev.off()
  
  ##################################################################
  # plot MSE
  ##################################################################
  
  mseVals<-c(res[,"cLRmse"],res[,"cLASSOmse"],res[,"pLRmse"],res[,"pLASSOmse"],res[,"p2LRmse"],res[,"p2LASSOmse"])
  
  pdf(paste("plotMSE","nA",nAdd,"nD",nDom,"nC",nCod,"nR",nRec,
            "MAFa",MAFadd*100,"MAFd",MAFdom*100,"MAFc",MAFcod*100,"MAFr",MAFrec*100,".pdf",sep=""))
  plot(-1,-1,xlim=c(min(betaS),max(betaS)),ylim=c(min(mseVals)-1*sd(mseVals),max(mseVals)+2*sd(mseVals)),ylab="Average MSE",
       xlab="Genetic Effect Size",cex=1.9)
  lines(betaS,res[,"cLRmse"],col=1,lty=3,type="b", pch=2,lwd=2)
  lines(betaS,res[,"cLASSOmse"],col=3,lty=3,type="b", pch=3,lwd=2)
  lines(betaS,res[,"pLRmse"],col=4,lty=3,type="b", pch=4,lwd=2)
  lines(betaS,res[,"pLASSOmse"],col=6,lty=3,type="b", pch=5,lwd=2)
  lines(betaS,res[,"p2LRmse"],col=7,lty=3,type="b", pch=6,lwd=2)
  lines(betaS,res[,"p2LASSOmse"],col=8,lty=3,type="b", pch=7,lwd=2)
  legend("topleft", inset=0.01,c("Classic PRS: regression","Classic PRS: LASSO",
                                 "Proposed PRS 1: regression","Proposed PRS 1: LASSO",
                                 "Proposed PRS 2: regression","Proposed PRS 2: LASSO"),col=c(1,3,4,6,7,8),pch=c(2:7),lty=3,lwd=2,cex=1)
  dev.off()
  
  ##################################################################
  # plot significant p-values
  ##################################################################
  
  pdf(paste("plotP","nA",nAdd,"nD",nDom,"nC",nCod,"nR",nRec,
            "MAFa",MAFadd*100,"MAFd",MAFdom*100,"MAFc",MAFcod*100,"MAFr",MAFrec*100,".pdf",sep=""))
  plot(-1,-1,xlim=c(min(betaS),max(betaS)),ylim=c(0,1),ylab="Proportion of simulations",
       xlab="Genetic Effect Size",cex=1.9)
  lines(betaS,res[,"cLRp"],col=1,lty=3,type="b", pch=2,lwd=2)
  lines(betaS,res[,"cLASSOp"],col=3,lty=3,type="b", pch=3,lwd=2)
  lines(betaS,res[,"pLRp"],col=4,lty=3,type="b", pch=4,lwd=2)
  lines(betaS,res[,"pLASSOp"],col=6,lty=3,type="b", pch=5,lwd=2)
  lines(betaS,res[,"p2LRp"],col=7,lty=3,type="b", pch=6,lwd=2)
  lines(betaS,res[,"p2LASSOp"],col=8,lty=3,type="b", pch=7,lwd=2)
  legend("topleft", inset=0.01,c("Classic PRS: regression","Classic PRS: LASSO",
                                 "Proposed PRS 1: regression","Proposed PRS 1: LASSO",
                                 "Proposed PRS 2: regression","Proposed PRS 2: LASSO"),col=c(1,3,4,6,7,8),pch=c(2:7),lty=3,lwd=2,cex=1)
  dev.off()
  }
    return(res)
}

