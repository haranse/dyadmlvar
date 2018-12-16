#' Analyze a Dyadic Timeseries Network.
#'
#' \code{read_network} Function that analyzes a dyadic timeseries network and calculates strength
#' and density variables. Assumes variables are list for one subject and then
#' equivalent list for another.
read_network <- function(mlvar_model) {
  ntwrk_vars <- list()
  ntwrk_vars$contemp <- list()
  ntwrk_vars$temporal <- list()


  #calculate node strength, as mean edge strength. Temporal node strength includes edges in both directions
  ntwrk_vars$contemp$strength <- rowMeans(abs(mlvar_model$results$Theta$cor$mean))-(1/length(mlvar_model$results$Theta$cor$mean[1,]))
  ntwrk_vars$temporal$strength <- rowMeans(abs(mlvar_model$results$Beta$mean))

  ntwrk_vars$vars <- list()
  ntwrk_vars$vars$names <- mlvar_model$input$vars
  ntwrk_vars$vars$len <- length(ntwrk_vars$vars$names)

  #getting couple-level indices
  ntwrk_vars$vars$num <-length(mlvar_model$results$Beta$subject[[1]])
  ntwrk_vars$temporal$all <- data.frame(matrix(nrow=0,ncol=ntwrk_vars$vars$num+1+sqrt(ntwrk_vars$vars$num)*2+3))
  ntwrk_vars$contemp$all <- data.frame(matrix(nrow=0,ncol=(ntwrk_vars$vars$num-sqrt(ntwrk_vars$vars$num))/2+3))

  for (i in 1:length(mlvar_model$IDs))
  #i = 1
     {
       temporal.1 <- mlvar_model$results$Beta$subject[[i]]
       dim(temporal.1)<-NULL
       temporal.1<-data.frame(temporal.1)
       temporal.1<-as.data.frame(t(temporal.1))

       #calculate strength either for connections to each node or from each node

       for (t in 1:(ntwrk_vars$vars$len)){
         temporal.1[[paste(ntwrk_vars$vars$names[[t]]," temp.to. centrality",sep="")]] <- rowMeans(abs(mlvar_model$results$Beta$subject[[i]]))[t]
         temporal.1[[paste(ntwrk_vars$vars$names[[t]]," temp.from. centrality",sep="")]] <- colMeans(abs(mlvar_model$results$Beta$subject[[i]]))[t]
       }

       #calculate per person density
       temporal.1.1<-data.frame(mlvar_model$results$Beta$subject[[i]])
       a<-mean(as.matrix((abs(temporal.1.1[1:(ntwrk_vars$vars$len/2),1:(ntwrk_vars$vars$len/2)]))))
       b<-mean(as.matrix((abs(temporal.1.1[(ntwrk_vars$vars$len/2+1):ntwrk_vars$vars$len,(ntwrk_vars$vars$len/2+1):ntwrk_vars$vars$len]))))
       c<-mean(as.matrix((abs(temporal.1.1[(ntwrk_vars$vars$len/2+1):ntwrk_vars$vars$len,1:(ntwrk_vars$vars$len/2)]))))
       d<-mean(as.matrix((abs(temporal.1.1[1:(ntwrk_vars$vars$len/2),(ntwrk_vars$vars$len/2+1):ntwrk_vars$vars$len]))))

       #intra density - connections in both directions between couple member's variables to themselves
       temporal.1$temporal.intra.density<-mean(c(a,b))

       #inter density - connections in both directions between couple member's variables to themselves
       temporal.1$temporal.inter.density<-mean(c(c,d))

       #ratio between intra and inter density
       temporal.1$temporal.ratio.density<-with(temporal.1,temporal.intra.density/temporal.inter.density)
       temporal.1$ID<-mlvar_model$IDs[i]
       ntwrk_vars$temporal$all<-rbind( ntwrk_vars$temporal$all,temporal.1)


       #Contemporaneous matrix is symmetrical - serializing half of it to avoid repetition of each edge twice
       contemporaneous.1 <- mlvar_model$results$Theta$pcor$subject[[i]]
       contemporaneous.2<-NULL
       for (t in 1:(ntwrk_vars$vars$len-1)){
         contemporaneous.2<-c(contemporaneous.2,contemporaneous.1[t,(t+1):ntwrk_vars$vars$len])
                                        }
       contemporaneous.2<-data.frame(t(contemporaneous.2))
       for (t in 1:(ntwrk_vars$vars$len)){
         contemporaneous.2[[paste(ntwrk_vars$vars$names[[t]]," cont. centrality",sep="")]] <- rowMeans(abs(mlvar_model$results$Theta$pcor$subject[[i]]))[t]-(1/length(mlvar_model$results$Theta$pcor$subject[[i]][1,]))
       }

       #caluculating inter and intra density
       cont.1<-data.frame(mlvar_model$results$Theta$pcor$subject[[i]])
       a<-sum(as.matrix((abs(cont.1[1:(ntwrk_vars$vars$len/2),1:(ntwrk_vars$vars$len/2)]))))
       a<-(a-(ntwrk_vars$vars$len/2))/((ntwrk_vars$vars$len/2)^2-(ntwrk_vars$vars$len/2))
       b<-sum(as.matrix((abs(cont.1[(ntwrk_vars$vars$len/2+1):ntwrk_vars$vars$len,(ntwrk_vars$vars$len/2+1):ntwrk_vars$vars$len]))))
       b<-(b-(ntwrk_vars$vars$len/2))/((ntwrk_vars$vars$len/2)^2-(ntwrk_vars$vars$len/2))
       c<-mean(as.matrix((abs(cont.1[(ntwrk_vars$vars$len/2+1):ntwrk_vars$vars$len,1:(ntwrk_vars$vars$len/2)]))))
       contemporaneous.2$cont.intra.density<-mean(c(a,b))
       contemporaneous.2$cont.inter.density<-c
       contemporaneous.2$cont.ratio.density<-with(contemporaneous.2,cont.intra.density/cont.inter.density)
       contemporaneous.2$ID<-mlvar_model$IDs[i]
       ntwrk_vars$contemp$all<-rbind(ntwrk_vars$contemp$all,contemporaneous.2)

  }

  #generate names for couple level variables
  ntwrk_vars$temporal$names <- vector()
  ntwrk_vars$temporal$centrality_names <- list()
  ntwrk_vars$temporal$inter_names <- vector()
  ntwrk_vars$temporal$intra_men_names <- vector()
  ntwrk_vars$temporal$intra_women_names <- vector()
  ntwrk_vars$temporal$intra_names <- vector()
  ntwrk_vars$temporal$density_names<-c("temporal.intra.density","temporal.inter.density","temporal.ratio.density")
  for(i in 1:ntwrk_vars$vars$len){
    for (j in 1:ntwrk_vars$vars$len)
    {
      ntwrk_vars$temporal$names[[(i-1)*ntwrk_vars$vars$len + j]]<-paste(ntwrk_vars$vars$names[[i]],"->",ntwrk_vars$vars$names[[j]])
    }
  }
  names(ntwrk_vars$temporal$all)[1:ntwrk_vars$vars$len^2] <- ntwrk_vars$temporal$names

  for (i in 1:length(ntwrk_vars$vars$names)){
    ntwrk_vars$temporal$centrality_names[[i*2 - 1]]<- paste(ntwrk_vars$vars$names[[i]]," temp.to. centrality",sep="")
    ntwrk_vars$temporal$centrality_names[[i*2]]<- paste(ntwrk_vars$vars$names[[i]]," temp.from. centrality",sep="")
  }
  half.len<-ntwrk_vars$vars$len/2
  for(i in 1:half.len){
    for (j in 1:half.len)
    {
      ntwrk_vars$temporal$intra_men_names[[(i-1)*half.len + j]]<-paste(ntwrk_vars$vars$names[[i]],"->",ntwrk_vars$vars$names[[j]])
      ntwrk_vars$temporal$inter_names[[(i-1)*half.len + j]]<-paste(ntwrk_vars$vars$names[[i]],"->",ntwrk_vars$vars$names[[j+half.len]])
    }

  }

  for(i in 1:half.len){
  for (j in 1:half.len)
  {
    ntwrk_vars$temporal$intra_women_names[[(i-1)*half.len + j]]<-paste(ntwrk_vars$vars$names[[i+half.len]],"->",ntwrk_vars$vars$names[[j+half.len]])
    ntwrk_vars$temporal$inter_names[[half.len^2 +(i-1)*half.len + j]]<-paste(ntwrk_vars$vars$names[[i+half.len]],"->",ntwrk_vars$vars$names[[j]])
  }
  }
  ntwrk_vars$temporal$intra_names<-c(ntwrk_vars$temporal$intra_men_names,ntwrk_vars$temporal$intra_women_names)

  ntwrk_vars$contemp$names <- vector()
  ntwrk_vars$contemp$centrality_names <- list()
  ntwrk_vars$contemp$inter_names <- vector()
  ntwrk_vars$contemp$intra_names <- vector()
  ntwrk_vars$contemp$intra_men_names <- vector()
  ntwrk_vars$contemp$intra_women_names <- vector()
  ntwrk_vars$contemp$density_names<-c("cont.intra.density","cont.inter.density","cont.ratio.density")
  for(i in 1:(ntwrk_vars$vars$len-1)){
    for (j in (i+1):ntwrk_vars$vars$len)
    {
      ntwrk_vars$contemp$names[(length(ntwrk_vars$contemp$names)+1)]<-paste(ntwrk_vars$vars$names[[i]],"<->",ntwrk_vars$vars$names[[j]])
    }
  }
  names(ntwrk_vars$contemp$all)[1:length(ntwrk_vars$contemp$names)]<-ntwrk_vars$contemp$names

  for (i in 1:length(ntwrk_vars$vars$names)){
    ntwrk_vars$contemp$centrality_names[i]<- paste(ntwrk_vars$vars$names[[i]]," cont. centrality",sep="")
  }

  for(i in 1:(half.len-1)){
    for (j in (i+1):half.len)
    {
      ntwrk_vars$contemp$intra_men_names[[length(ntwrk_vars$contemp$intra_men_names)+1]]<-paste(ntwrk_vars$vars$names[[i]],"<->",ntwrk_vars$vars$names[[j]])
      ntwrk_vars$contemp$intra_women_names[[length(ntwrk_vars$contemp$intra_women_names)+1]]<-paste(ntwrk_vars$vars$names[[i+half.len]],"<->",ntwrk_vars$vars$names[[j+half.len]])
          }

  }
  ntwrk_vars$contemp$intra_names<-c(ntwrk_vars$contemp$intra_men_names,ntwrk_vars$contemp$intra_women_names)
  for(i in 1:(half.len)){
    for (j in (half.len+1):(half.len*2))
    {
      ntwrk_vars$contemp$inter_names[[length(ntwrk_vars$contemp$inter_names)+1]]<-paste(ntwrk_vars$vars$names[[i]],"<->",ntwrk_vars$vars$names[[j]])

    }


  }

  return(ntwrk_vars)
}

lasso<-function(Data,Predictors,Outcome,Seeds=1,Train=F,PropOfTrain=.75){
  #Data=all,Predictors=pre_inter,Outcome="W_csi_3_resid_1",Train=F,PropOfTrain=.75)
  all.2<-Data[,c(Predictors,Outcome)]
  #removing incomplete data;
  all.2<-stats::na.omit(all.2)
  #standardized data- #see here for more info: https://stats.stackexchange.com/questions/126109/coefficient-value-from-glmnet
  #and here https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html
  all.2<-as.data.frame(scale(all.2,center=T,scale = T))
  #taking out edges without variability
  for (no.var in Predictors){
    #no.var<-Predictors[9]
    print(paste(no.var,":",is.na(stats::var(all.2[,no.var]))))
    print(stats::var(all.2[,no.var]))
    if (is.na(stats::var(all.2[,no.var]))){
      all.2<-all.2[,-which(names(all.2)==no.var)]
        }
  }
  #transforming into a matrix
  x<-stats::model.matrix (all.2[,Outcome]~.,all.2)[,-1]
  x<-x[,-(length(x[1,]))]
  y<-as.matrix(all.2[,Outcome])
  #If  train=F then no splitting of data occurs, otherwise data is splitted
  if (Train==F){
    #using cross-validation to choose the tuning parameter; by default it performs ten-fold cross-validation
    set.seed(Seeds)
    cv.out<-glmnet::cv.glmnet(x,y,alpha=1,intercept = F,standardize = F)
    graphics::plot(cv.out)
    #for plot interpretatin see: https://stats.stackexchange.com/questions/253963/how-to-interpret-cv-glmnet-plot
    bestlam <-cv.out$lambda.min
    #getting estimates
    a<-glmnet::glmnet(x,y,alpha=1,lambda = bestlam,intercept = F,standardize = F)
    b<-stats::coefficients(a)
    b<-b[1:length(b),]
    b<-(round(b[b!=0],3))
    #estimate a ranking of variable importance,
    #i.e., the maximal value of lambda at which
    #the variable first entered the model was determined
    d<-as.data.frame(b)
    names(d)[1]<-"Est."
    d.2<-rownames(d)
    for (i in d.2){
      #i<-d.2[1]
      e<-as.data.frame(cv.out$glmnet.fit$beta[i,])
      row.names(e)<-NULL
      e<-subset(e,e[,1]!=0)
      d[i,"Importance"]<-row.names(e)[1]
    }

    pred.mat <- cbind(y_true = y, y_pred = stats::predict(a,newx=x)[, 1])
    Rs <- boot::corr(d=pred.mat) ^ 2
    c<-list(Estimates=d,Rsquared=Rs)
    return(c)
  } else {
        #splitting into train and test data
  set.seed(Seeds)
  train<-sort(sample (1: nrow(x), nrow(x)*PropOfTrain))
  test<-(-train)
  y.test<-y[test]
  #using cross-validation to choose the tuning parameter; by default it performs ten-fold cross-validation
  #At least 10 observation per fold
  folds<-floor(length(y[train])/10)
  set.seed(Seeds)
  cv.out<-glmnet::cv.glmnet(x[train,],y[train],alpha=1,intercept = F,standardize = F,nfolds = folds)
  graphics::plot(cv.out)
  #for plot interpretatin see: https://stats.stackexchange.com/questions/253963/how-to-interpret-cv-glmnet-plot
  bestlam <-cv.out$lambda.min
  #getting estimates
  a<-glmnet::glmnet(x[train,],y[train],alpha=1,lambda = bestlam,intercept = F,standardize = F)
  b<-stats::coefficients(a)
  b<-b[1:length(b),]
  b<-(round(b[b!=0],3))
  #estimate a ranking of variable importance,
  #i.e., the maximal value of lambda at which
  #the variable first entered the model was determined
  d<-as.data.frame(b)
  names(d)[1]<-"Est."
  d.2<-rownames(d)
    for (i in d.2){
    #i<-d.2[1]
    e<-as.data.frame(cv.out$glmnet.fit$beta[i,])
    row.names(e)<-NULL
    e<-subset(e,e[,1]!=0)
    d[i,"Importance"]<-row.names(e)[1]
  }

  pred.mat <- cbind(y_true = y[test], y_pred = stats::predict(a,newx=x[test,])[, 1])
  Rs <- boot::corr(d=pred.mat) ^ 2
  c<-list(Estimates=d,Rsquared=Rs)
  return(c)
  }
}
