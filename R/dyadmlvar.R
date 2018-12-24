#' @details The package analyzes dyadic networks computed by the mlVAR package.
#'   Use read_network to analyze a network, get_names to get lists of variable
#'   names for specific types of variables (e.g. inter-partner variables), and
#'   subset them from the table of all variables (ntwrk$all if ntwrk is the
#'   output of read_network). Due to the large number of variables you will
#'   probably want to use the lasso function to see which variables explain
#'   significant amounts of variance in some outcome vector. More detail can be
#'   found in Bar-Kalifa & Sened, 2019 (Under review).
"_PACKAGE"

#' Example mlVAR model
#' @description An mlVAR model obtained from diary data collected from
#'   80 couples. Nodes include daily feelings of anger, sadness and anxiety for
#'   each partner an mlVAR object, with some internal variables removed to
#'   conserve package space
#' @format an mlVAR object with 6 nodes
#' @source an internal dataset, analyzed using the mlVAR package (Epskamp,
#'   Deserno & Bringmann, 2018). For details on the data collection see:
#'   Bar-Kalifa, E., Rafaeli, E., & Sened, H. (2016). Truth and bias in daily
#'   judgments of support receipt between romantic partners. Personal
#'   Relationships, 23(1), 42-61.
"fit1"

#' Relationship satisfaction data
#' @description relationship satisfaction before and after completing a daily
#'   diary for 80 couples, and residuals of post-diary satisfaction after
#'   adjusting for pre-diary satisfaction.
#' @format a dataframe with 7 variables and 80 rows
#' @source an internal dataset, for details on the data collection see:
#'   Bar-Kalifa, E., Rafaeli, E., & Sened, H. (2016). Truth and bias in daily
#'   judgments of support receipt between romantic partners. Personal
#'   Relationships, 23(1), 42-61.
"sat"


#' Analyze a Dyadic Timeseries Network.
#'
#' \code{read_network} Function that analyzes a dyadic timeseries network and
#' calculates strength and density variables. Assumes variables are list for one
#' subject and then equivalent list for another.
#' @param mlvar_model An mlVAR object, usually created by the mlVAR function of
#'   package mlVAR. The target variables should have been of even length, with
#'   the first half belonging to part A of the dyad (e.g. men, clients,
#'   children) and the other to part B (e.g. women, therapists, parents)
#' @return a dyadNetwork object with variables characterizing the network
#' @export
#' @examples
#' ntwrk <- read_network(fit1)
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
         temporal.1[[paste(ntwrk_vars$vars$names[[t]]," temp.to. strength",sep="")]] <- rowMeans(abs(mlvar_model$results$Beta$subject[[i]]))[t]
         temporal.1[[paste(ntwrk_vars$vars$names[[t]]," temp.from. strength",sep="")]] <- colMeans(abs(mlvar_model$results$Beta$subject[[i]]))[t]
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
         contemporaneous.2[[paste(ntwrk_vars$vars$names[[t]]," cont. strength",sep="")]] <- rowMeans(abs(mlvar_model$results$Theta$pcor$subject[[i]]))[t]-(1/length(mlvar_model$results$Theta$pcor$subject[[i]][1,]))
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
  ntwrk_vars$temporal$strength_names <- list()
  ntwrk_vars$temporal$inter_names <- vector()
  ntwrk_vars$temporal$intra_partA_names <- vector()
  ntwrk_vars$temporal$intra_partB_names <- vector()
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
    ntwrk_vars$temporal$strength_names[[i*2 - 1]]<- paste(ntwrk_vars$vars$names[[i]]," temp.to. strength",sep="")
    ntwrk_vars$temporal$strength_names[[i*2]]<- paste(ntwrk_vars$vars$names[[i]]," temp.from. strength",sep="")
  }
  half.len<-ntwrk_vars$vars$len/2
  for(i in 1:half.len){
    for (j in 1:half.len)
    {
      ntwrk_vars$temporal$intra_partA_names[[(i-1)*half.len + j]]<-paste(ntwrk_vars$vars$names[[i]],"->",ntwrk_vars$vars$names[[j]])
      ntwrk_vars$temporal$inter_names[[(i-1)*half.len + j]]<-paste(ntwrk_vars$vars$names[[i]],"->",ntwrk_vars$vars$names[[j+half.len]])
    }

  }

  for(i in 1:half.len){
  for (j in 1:half.len)
  {
    ntwrk_vars$temporal$intra_partB_names[[(i-1)*half.len + j]]<-paste(ntwrk_vars$vars$names[[i+half.len]],"->",ntwrk_vars$vars$names[[j+half.len]])
    ntwrk_vars$temporal$inter_names[[half.len^2 +(i-1)*half.len + j]]<-paste(ntwrk_vars$vars$names[[i+half.len]],"->",ntwrk_vars$vars$names[[j]])
  }
  }
  ntwrk_vars$temporal$intra_names<-c(ntwrk_vars$temporal$intra_partA_names,ntwrk_vars$temporal$intra_partB_names)

  ntwrk_vars$contemp$names <- vector()
  ntwrk_vars$contemp$strength_names <- list()
  ntwrk_vars$contemp$inter_names <- vector()
  ntwrk_vars$contemp$intra_names <- vector()
  ntwrk_vars$contemp$intra_partA_names <- vector()
  ntwrk_vars$contemp$intra_partB_names <- vector()
  ntwrk_vars$contemp$density_names<-c("cont.intra.density","cont.inter.density","cont.ratio.density")
  for(i in 1:(ntwrk_vars$vars$len-1)){
    for (j in (i+1):ntwrk_vars$vars$len)
    {
      ntwrk_vars$contemp$names[(length(ntwrk_vars$contemp$names)+1)]<-paste(ntwrk_vars$vars$names[[i]],"<->",ntwrk_vars$vars$names[[j]])
    }
  }
  names(ntwrk_vars$contemp$all)[1:length(ntwrk_vars$contemp$names)]<-ntwrk_vars$contemp$names

  for (i in 1:length(ntwrk_vars$vars$names)){
    ntwrk_vars$contemp$strength_names[i]<- paste(ntwrk_vars$vars$names[[i]]," cont. strength",sep="")
  }

  for(i in 1:(half.len-1)){
    for (j in (i+1):half.len)
    {
      ntwrk_vars$contemp$intra_partA_names[[length(ntwrk_vars$contemp$intra_partA_names)+1]]<-paste(ntwrk_vars$vars$names[[i]],"<->",ntwrk_vars$vars$names[[j]])
      ntwrk_vars$contemp$intra_partB_names[[length(ntwrk_vars$contemp$intra_partB_names)+1]]<-paste(ntwrk_vars$vars$names[[i+half.len]],"<->",ntwrk_vars$vars$names[[j+half.len]])
          }

  }
  ntwrk_vars$contemp$intra_names<-c(ntwrk_vars$contemp$intra_partA_names,ntwrk_vars$contemp$intra_partB_names)
  for(i in 1:(half.len)){
    for (j in (half.len+1):(half.len*2))
    {
      ntwrk_vars$contemp$inter_names[[length(ntwrk_vars$contemp$inter_names)+1]]<-paste(ntwrk_vars$vars$names[[i]],"<->",ntwrk_vars$vars$names[[j]])

    }


  }
  ntwrk_vars$all <- base::merge(ntwrk_vars$temporal$all, ntwrk_vars$contemp$all, by="ID")
  class(ntwrk_vars) <- "dyadNetwork"
  return(ntwrk_vars)
}

#' Constant for retreiving node connections between partner A variables
#' @export
INTRA_A = 1

#' Constant for retreiving node connections between partner B variables
#' @export
INTRA_B = 2

#' Constant for retreiving node connections between partner A variables to
#' parter B variables
#' @export
INTER = 3

#' Constant for retreiving total node strengths
#' @export
STRENGTH = 4

#' Constant for retreiving densities of intra-partner networks, inter-partner
#' networks, and their ratio
#' @export
DENSITY = 5

#' Constant for retreiving variables from the contemporaneous network
#' @export
CONTEMP = 1

#' Constant for retreiving variables from the temporal network
#' @export
TEMPORAL = 2

#' Constant for retreiving all available variables
#' @export
ALL_NETWORK = 6

#' Get a list of variable names form a dyadNetwork object
#'
#' \code{get_names} Get a list of variable names form a dyadNetwork object,
#' e.g. variable names for intra-partner connections for partner A
#' @param ntwrk a dyadNetwork object
#' @param part which kind of variables to retrieve:
#'      INTRA_A - interconnections between partner A variables
#'      INTRA_B - interconnections between partner B variables
#'      INTER - connections between partner A to partner B variables
#'      STRENGTH - node strengths
#'      DENSITY - density of intra-partner and inter-partner networks and
#'          the ratio between both densities
#'      ALL_NETWORK - all variables
#' @param time which network to retreive variables from:
#'      TEMPORAL - temporal network
#'      CONTEMP - contemporaneous network
#'      ALL_NETWORK - both
#' @return a list of variable name strings
#' @export
#' @examples
#' ntwrk <- read_network(fit1)
#' #all variable names
#' all_vars <- get_names(ntwrk)
#' #inter-partner variables from the contemporaneous network
#' inter_contemp_vars <- get_names(ntwrk, part = INTER, time = CONTEMP)
get_names<-function(ntwrk, part = ALL_NETWORK, time = ALL_NETWORK) {
  if (time == ALL_NETWORK) {
    return(c(get_names(ntwrk,part,CONTEMP),get_names(ntwrk,part,TEMPORAL)))
  }
  if (time == CONTEMP)
  {
    target = ntwrk$contemp
  }
  else
  {
    target = ntwrk$temporal
  }

  if (part == ALL_NETWORK)
  {
    return(c(get_names(ntwrk,INTRA_A,time),get_names(ntwrk,INTRA_B,time),get_names(ntwrk,INTER,time),
             get_names(ntwrk,STRENGTH,time),get_names(ntwrk,DENSITY,time)))
  }
  else if (part == INTRA_A)
  {
    return(target$intra_partA_names)
  }
  else if (part == INTRA_B)
  {
    return(target$intra_partB_names)
  }
  else if (part == INTER)
  {
    return(target$inter_names)
  }
  else if (part == STRENGTH)
  {
    return(target$strength_names)
  }
  else if (part == DENSITY)
  {
     return(target$density_names)
  }
}


#' Use a LASSO algorithm to find important predictors
#'
#' \code{lasso} Function that uses the LASSO shrinkage reduction method
#' to find meaningful predictors of an outcome vector
#' @param Data a dataframe with columns for each predictor and for the outcome
#'   variable
#' @param Predictors a list of strings with predictor variables names
#' @param Outcome a string with the name of an outcome variable
#' @param Seeds a seed for randomly selecting some part of the data to be
#'   training data and some to be testing data. default -
#' @param Train Should we split data to training and test datasets (FALSE uses
#'   all data for both)
#' @param PropOfTrain How much of the data to use for training the model
#' @param Plot shoud we show a plot of the parameter number and log likelihood?
#' @return a list of estimates for the effect of the valid predictors, and a R^2
#'   statistic for the final model
#' @export
#' @examples
#' ntwrk <- read_network(fit1)
#' #intra-partner variables from the partner B
#' intra_B_vars <- get_names(ntwrk, part = INTRA_B, time = ALL_NETWORK)
#' full_data <- merge(ntwrk$all, sat, by = "ID", all=TRUE, sort=TRUE)
#' lasso(full_data, intra_B_vars, "W_csi_resid")
lasso<-function(Data,Predictors,Outcome,Seeds=1,Train=F,PropOfTrain=.75, Plot=F){
  all.2<-Data[,c(Predictors,Outcome)]
  #removing incomplete data;
  all.2<-stats::na.omit(all.2)
  #standardized data- #see here for more info:
  #https://stats.stackexchange.com/questions/126109/coefficient-value-from-glmnet
  #and here https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html
  all.2<-as.data.frame(scale(all.2,center=T,scale = T))
  #taking out edges without variability
  for (no.var in Predictors){
    #no.var<-Predictors[9]
    #print(paste(no.var,":",is.na(stats::var(all.2[,no.var]))))
    #print(stats::var(all.2[,no.var]))
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
    if (Plot)
    {
      graphics::plot(cv.out)
    }
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
  if (Plot)
  {
    graphics::plot(cv.out)
  }
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
