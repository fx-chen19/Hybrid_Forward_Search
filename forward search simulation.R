library(ts.extend)
library(openxlsx)
#Functions-------------------
#Create dataset
set_ds <- function(N,T,n_clusters
                   ,with_epidemic,change_param,start_tpt,affected_period,n_clusters_affected){
  t <- sort(rep(1:N,T))
  space <- rep(1:N,T)
  n = N / n_clusters
  
  w1 <- rpois(T*n,lambda_c1)
  w2 <- rpois(T*n,lambda_c2)
  w3 <- rpois(T*n,lambda_c3)
  w4 <- rpois(T*n,lambda_c4)
  w5 <- rpois(T*n,lambda_c5)

  cluster <- c(); w <- c()
  for (i in 1:n_clusters){
    cluster <- c(cluster,rep(i,floor(n)))
    nam <- paste0("w", i)
    assign(nam,get(nam))
    w <- append(w,get(nam))
  }
  
  #handling for unbalanced dataset
  len=length(cluster)
  if (len!=N){
    for (i in c(1:(N-len))){
      cluster <- append(cluster,i)
    }
  }
  cluster <- sort(cluster)
  cluster <- rep(cluster,T)
  
  if (with_epidemic==FALSE){
    beta=beta; gamma=gamma; lambda0=0
    surge=0
  } else{ #with structural change
    affected_clusters <- sample(c(1:n_clusters), n_clusters_affected,replace=FALSE)
    affected_tpt = floor(affected_period*T)
    surge <- ifelse(t<(start_tpt+affected_tpt) & t>=start_tpt & cluster %in% affected_clusters,1,0)
    beta_s=beta*(1+change_param)
    gamma_s=gamma*(1+change_param)
    gamma <- ifelse(t<(start_tpt+affected_tpt) & t>=start_tpt & cluster %in% affected_clusters,gamma_s,gamma)
    beta <- ifelse(t<(start_tpt+affected_tpt) & t>=start_tpt & cluster %in% affected_clusters,beta_s,beta)
    lambda0 <- ifelse(t<(start_tpt+affected_tpt) & t>=start_tpt & cluster %in% affected_clusters,lambda0,0)
  }
  # X~N(mean=10,000, var = 1,000) = non neighborhood vars
  x <- rnorm(N*T, mean=10000, sd=sqrt(1000))
  
  #AR(1)
  at <- rGARMA(n=N, m=T, ar=rho, mean=0, errorvar=1)
  at <- matrix(at,nrow=N*T,byrow=FALSE)
  
  y <- beta*x+gamma*w+at+lambda0*exp(-lambda1*t)
  
  df <- cbind(t,space,cluster,surge,y,x,w)
  colnames(df)[5]='y'
  df <- as.data.frame(df)
  
  print('df created')
  return(df)
}

#get rho and mad
rho_mad <- function(fwddata2){
  #step2 nonepi#
  #AR on residuals to get rho
  
  ARfwd<-fwddata2
  dum<-ARfwd$resid
  dum<-c(NA,dum)
  dum<-head(dum,-1) #remove last observation
  ARfwd[ARfwd$t==1,'resid'] <-NA #resid col

    #get average for each area
  acflist<-data.frame(matrix(ncol=1,nrow=0))
  colnames(acflist)<-c('acf')
  
  for (k in unique(fwddata2$space)) {
    subprov <- subset.data.frame(fwddata2,space==k)
    acfest <- acf(subprov$resid, lag.max = 1, type = "correlation", plot = F)
    acflist<-rbind(acflist,acfest$acf[2])
  }

  colnames(acflist)<-c('rho')
  rho <- colMeans(acflist)

  #MAD
  ARfwd['MAD'] <- abs(ARfwd$y-ARfwd$yhat)
  MAD <- mean(ARfwd$MAD)

    #MAPE
  ARfwd['MAPE'] <- abs((ARfwd$y-ARfwd$yhat)/ARfwd$y)
  MAPE <- mean(ARfwd$MAPE)
    
  res <- data.frame(rho=rho,MAD=MAD,MAPE=MAPE)
  return(res)
}

#forward search construction
hybrid <- function(fwddata,N,T,n_clusters,with_epidemic){
  #step1
  fwd_lm1<-lm(y~x+w,data=fwddata)
  fwd_coeff=as.data.frame(t(fwd_lm1$coefficients))[,]
  
  for (i in unique(fwddata$t)) {
    
    subdata<-subset.data.frame(fwddata,t==i) #select data from each t point
    fwd_lm1<-lm(y~x+w,data=subdata)
    summary(fwd_lm1)
    
    residual<-fwd_lm1$residuals
    subdata<-cbind(subdata,residual)
    subdata<-subdata[order(subdata$residual),]
    subset1<-subdata[1:ceiling(nrow(subdata)/4),] #25%
    subset2<-subdata[ceiling(nrow(subdata)/4)+1:nrow(subdata),]  
    subset2<-subset2[complete.cases(subset2),] #75%
    ncook <- 1
    
    while (ncook > 0) {
      i=1
      fwd_lm2<-lm(y~x+w,data=subset1)
      if (nrow(subset1)==nrow(subdata)) {
        ncook <- 0
      }else {
        subset2$fitted.value<-predict(object = fwd_lm2,newdata = subset2)
        subset2$residual<-subset2$y-subset2$fitted.value
        
        subset2<-subset2[order(subset2$residual),] #order by resids
        addrow<-subset(subset2,select=-c(fitted.value)) #get obs with smallest residual and remove resids col
        
        subset1<-rbind(subset1,addrow)
        subset2<-subset2[-1,]
        
        cooksd <- cooks.distance(fwd_lm2)>4/nrow(subset1) #find rows in ols where cooks d is greater than 4/nrow(subset1)
        ncook <- sum(cooksd, na.rm = TRUE)
        i=i+1
        }
    }
    #The algorithm then stops if the Cook’s D is no longer influential to the model based on this threshold
    
    fwd_coeff<-rbind(fwd_coeff,fwd_lm2$coefficients)
  }
  fwd_coeff<-fwd_coeff[-1,] #remove test row
  fws_mncoeff<-as.data.frame(t(colMeans(fwd_coeff))) 
  
  fwddata2 <- fwddata
  fwddata2['yhat']<-  fws_mncoeff[1,'(Intercept)']+fws_mncoeff[1,'x']*fwddata['x']+fws_mncoeff[1,'w']*fwddata['w']
  fwddata2['resid']<-fwddata2$y-fwddata2$yhat
  
  ols <- fwddata
  ols_lm<-lm(y~x+w,data=fwddata)
  ols_coeff <- ols_lm$coefficients
  ols_coeff <- t(as.data.frame(ols_coeff))
  ols['yhat']<- predict(object = ols_lm,newdata = ols)
  ols['resid']<-ols$y-ols$yhat
  
  if (with_epidemic==FALSE){
    stats <- rho_mad(fwddata2)
    olstats <- rho_mad(ols)
  } else{
    #lambda0 and lambda1 with epidemic
    fwddata.epi<-fwddata2
    fwdepi.surge<-subset.data.frame(fwddata.epi,surge==1)
    expo<-lm(log(resid)~t,data=fwdepi.surge) #may negative na residuals, no log
    print(expo)
    lambda0<-exp(expo$coefficient[1])
    lambda1<- -expo$coefficient[2]
    
    fwddata.epi['yhat']<-
      fws_mncoeff[1,'(Intercept)']+fws_mncoeff[1,'x']*fwddata['x']+fws_mncoeff[1,'w']*fwddata['w']+
      fwddata$surge*lambda0*exp(-lambda1*fwddata$t)
    fwddata.epi['resid'] <-fwddata.epi$y-fwddata.epi$yhat
    print (fwddata.epi)
    
    
    print('fit and resids')
    epi_stats <- rho_mad(fwddata.epi)
    print('rho obtained')
      #OLS
    ols.surge<-subset.data.frame(ols,surge==1)
    ols.surge
    expo<-lm(log(resid)~t,data=ols.surge) #may negative na residuals, no log
    print(expo)
    lambda0_ols<-exp(expo$coefficient[1])
    lambda1_ols<- -expo$coefficient[2]
    ols['yhat']<-  
      ols_coeff[1,'(Intercept)']+ols_coeff[1,'x']*fwddata['x']+ols_coeff[1,'w']*fwddata['w']+
          fwddata$surge*lambda0_ols*exp(-lambda1_ols*fwddata$t)
    ols['resid']<-ols$y-ols$yhat
    olstats <- rho_mad(ols)
  }    
  
  config <- data.frame(T=T,N=N,Clusters=n_clusters)
  if (with_epidemic==FALSE){
      est <- data.frame(beta=fws_mncoeff$x,beta_ols=ols_coeff[1,'x']
             ,gamma=fws_mncoeff$w,gamma_ols=ols_coeff[1,'w']
             ,rho=stats$rho,rho_ols=olstats['rho']
             ,MAD=stats$MAD,MAD_ols=olstats['MAD']
             )
  }else{
      est <- data.frame(beta=fws_mncoeff$x,beta_ols=ols_coeff[1,'x']
             ,gamma=fws_mncoeff$w,gamma_ols=ols_coeff[1,'w']
             ,lambda0=lambda0,lambda0_ols=lambda0_ols
             ,lambda1=lambda1,lambda1_ols=lambda1_ols
             ,rho=epi_stats$rho,rho_ols=olstats['rho']
             ,MAD=epi_stats$MAD,MAD_ols=olstats['MAD']
             )
  }    
  
  fwd_finalcoeff <- cbind(config,est)
  return(fwd_finalcoeff)
}
execute_program <- function(T,N,n_clusters=c(2,5),with_epidemic
                            ,change_param=0,start_tpt=0,affected_period=0
                            ,n_clusters_affected=0
                            ,case=''
                            ,sheetname){
  final <- data.frame()
  for (change in change_param){
  for (ti in T){
    for (ni in N){
      for (n_cluster in n_clusters){ #repeat for each cluster configuration
        if (tolower(n_clusters_affected)=='all'){
          n_clusters_affected=n_cluster
        }
        df <- set_ds(N=ni,T=ti,n_clusters=n_cluster
                     ,with_epidemic,change_param=change,start_tpt,affected_period
                     ,n_clusters_affected
                     )
        hybrid_df <- hybrid(df,ni,ti,n_cluster,with_epidemic)
        final <- rbind(final, hybrid_df)
      }
    }
  }
  }
  names(final)[names(final) == 'rho.1'] <- 'rho_ols'
  names(final)[names(final) == 'MAD.1'] <- 'MAD_ols'
  if (with_epidemic){
    final <- cbind(Case=case,'Change in Param'=change_param,final)
  }
  
  wb <- loadWorkbook(file = "Estimates.xlsx")
  if (sheetname %in% sheets(wb)){
    writeData(wb, sheetname, final
              ,startRow = nrow(readWorkbook("Estimates.xlsx",sheet=sheetname))+2
              ,colNames=FALSE)
    
  } else{
    addWorksheet(wb, sheetname)
    writeData(wb, sheetname, final
              ,startRow = 1
              ,colNames=TRUE)
  }
  saveWorkbook(wb, file = "Estimates.xlsx", overwrite = TRUE)
  
  
  #MAPE
  rows=dim(final)[1]
  betas=rep(beta,rows)
  gammas=rep(gamma,rows)
  rhos=rep(rho,rows)
  lambda0s=rep(lambda0,rows)
  lambda1s=rep(lambda1,rows)
  
  if(with_epidemic==TRUE){
  percent_diff <- final[1:4]
  percent_diff$beta <- abs((final$beta-betas)/betas)
  percent_diff$beta_ols <- abs((final$beta_ols-betas)/betas)
  percent_diff$gamma <- abs((final$gamma-gammas)/gammas)
  percent_diff$gamma_ols <- abs((final$gamma_ols-gammas)/gammas)
  percent_diff$lambda0 <- abs((final$lambda0-lambda0)/lambda0)
  percent_diff$lambda0_ols <- abs((final$lambda0_ols-lambda0)/lambda0)
  percent_diff$lambda1 <- abs((final$lambda1-lambda1)/lambda1)
  percent_diff$lambda1_ols <- abs((final$lambda1_ols-lambda1)/lambda1)
  percent_diff$rho <- abs((final$rho-rhos)/rhos)
  percent_diff$rho_ols <- abs((final$rho_ols-rhos)/rhos)
  mapes <- data.frame(MAPE=rowMeans(percent_diff[,c(5,7,9,11,13)]),MAPE_ols=rowMeans(percent_diff[,c(6,8,10,12,14)]))
  } else{
  percent_diff <- final[1:3]
  percent_diff$beta <- abs((final$beta-betas)/betas)
  percent_diff$beta_ols <- abs((final$beta_ols-betas)/betas)
  percent_diff$gamma <- abs((final$gamma-gammas)/gammas)
  percent_diff$gamma_ols <- abs((final$gamma_ols-gammas)/gammas)
  percent_diff$rho <- abs((final$rho-rhos)/rhos)
  percent_diff$rho_ols <- abs((final$rho_ols-rhos)/rhos)
  mapes <- data.frame(MAPE=rowMeans(percent_diff[,c(4,6,8)]),MAPE_ols=rowMeans(percent_diff[,c(5,7,9)]))
  }
  
  
  percent_diff <- cbind(percent_diff,mapes)

    if (with_epidemic){
      percent_diff <- cbind(Case=case,'Change in Param'=change_param,percent_diff)
  }
  
  wb <- loadWorkbook(file = "Percent Differences.xlsx")
  if (sheetname %in% sheets(wb)){
    writeData(wb, sheetname, percent_diff
              ,startRow = nrow(readWorkbook("Percent Differences.xlsx",sheet=sheetname))+2
              ,colNames=FALSE)
    
  } else{
    addWorksheet(wb, sheetname)
    writeData(wb, sheetname, percent_diff
              ,startRow = 1
              ,colNames=TRUE)
  }
  saveWorkbook(wb, file = "Percent Differences.xlsx", overwrite = TRUE)
}


#Program Execution===============
#set seed here to get consistent results
set.seed(0) 

#base parameters
beta=0.52; gamma=14.6; rho=0.5 
lambda0=4.8*10^6; lambda1=2.5
#w
lambda_c1=100 ;lambda_c2=200 ;lambda_c3=300 
lambda_c4=400 ;lambda_c5=500 

#create empty xlsx for data
write.xlsx('','Estimates.xlsx')
write.xlsx('','Percent Differences.xlsx')

## without epidemic----------------
#put diff configuration here
execute_program(T = c(10,20,50)
                ,N = c(20,30,50)
                ,with_epidemic = FALSE
                ,sheetname = 'non-epidemic case'
                )

execute_program(T = c(10,20)
                ,N = c(25,26,27)
                ,with_epidemic = FALSE
                ,sheetname = 'non-epidemic case'
)
##with epidemic---------------
### 1 contaminated-----------------
#### Case 1: contamination in 1 cluster, short period, no change in parameters ------
execute_program(T = 20
                ,N = 20
                ,with_epidemic = TRUE
                ,change_param = 0
                ,start_tpt = 1
                ,affected_period = 0.25
                ,n_clusters_affected = 1
                ,case='Case 1'
                ,sheetname = '1 Contaminated Small Balanced'
                )
execute_program(T = 50
                ,N = 50
                ,with_epidemic = TRUE
                ,change_param = 0
                ,start_tpt = 1
                ,affected_period = 0.25
                ,n_clusters_affected = 1
                ,case='Case 1'
                ,sheetname = '1 Contaminated Large Balanced'
                )
execute_program(T = 10
                ,N = 27
                ,with_epidemic = TRUE
                ,change_param = 0
                ,start_tpt = 1
                ,affected_period = 0.25
                ,n_clusters=2
                ,n_clusters_affected = 1
                ,case='Case 1'
                ,sheetname = '1 Contaminated Small Unbalanced'
)
execute_program(T = 10
                ,N = 27
                ,with_epidemic = TRUE
                ,change_param = 0
                ,start_tpt = 1
                ,affected_period = 0.25
                ,n_clusters=5
                ,n_clusters_affected = 1
                ,case='Case 1'
                ,sheetname = '1 Contaminated Small Unbalanced'
)
#### Case 2: contamination in 1 cluster, short period, with change in parameters ------
execute_program(T = 20
                ,N = 20
                ,with_epidemic = TRUE
                ,change_param = c(0.10,0.20,0.30,0.40)
                ,start_tpt = 1
                ,affected_period = 0.25
                ,n_clusters_affected = 1
                ,case='Case 2'
                ,sheetname = '1 Contaminated Small Balanced'
                )
execute_program(T = 50
                ,N = 50
                ,with_epidemic = TRUE
                ,change_param = c(0.10,0.20,0.30,0.40)
                ,start_tpt = 1
                ,affected_period = 0.25
                ,n_clusters_affected = 1
                ,case='Case 2'
                ,sheetname = '1 Contaminated Large Balanced'
                )
execute_program(T = 10
                ,N = 29
                ,with_epidemic = TRUE
                ,change_param = c(0.10,0.20,0.30,0.40)
                ,start_tpt = 1
                ,affected_period = 0.25
                ,n_clusters=2
                ,n_clusters_affected = 1
                ,case='Case 2'
                ,sheetname = '1 Contaminated Small Unbalanced'
                )
execute_program(T = 10
                ,N = 22
                ,with_epidemic = TRUE
                ,change_param = c(0.10,0.20,0.30,0.40)
                ,start_tpt = 1
                ,affected_period = 0.25
                ,n_clusters=5
                ,n_clusters_affected = 1
                ,case='Case 2'
                ,sheetname = '1 Contaminated Small Unbalanced'
                )


#### Case 3: contamination in 1 cluster, long period, no change in parameters ------
execute_program(T = 20
                ,N = 20
                ,with_epidemic = TRUE
                ,change_param = 0
                ,start_tpt = 1
                ,affected_period = 0.50
                ,n_clusters_affected = 1
                ,case='Case 3'
                ,sheetname = '1 Contaminated Small Balanced'
                )
execute_program(T = 10
                ,N = 27
                ,with_epidemic = TRUE
                ,change_param = 0
                ,start_tpt = 1
                ,affected_period = 0.50
                ,n_clusters=2
                ,n_clusters_affected = 1
                ,case='Case 3'
                ,sheetname = '1 Contaminated Small Unbalanced'
                )
execute_program(T = 10
                ,N = 29
                ,with_epidemic = TRUE
                ,change_param = 0
                ,start_tpt = 1
                ,affected_period = 0.50
                ,n_clusters=5
                ,n_clusters_affected = 1
                ,case='Case 3'
                ,sheetname = '1 Contaminated Small Unbalanced'
                )
execute_program(T = 50
                ,N = 50
                ,with_epidemic = TRUE
                ,change_param = 0
                ,start_tpt = 1
                ,affected_period = 0.50
                ,n_clusters_affected = 1
                ,case='Case 3'
                ,sheetname = '1 Contaminated Large Balanced'
                )
#### Case 4: contamination in 1 cluster, long period, with change in parameters ------
execute_program(T = 20
                ,N = 20
                ,with_epidemic = TRUE
                ,change_param = c(0.10,0.20,0.30,0.40)
                ,start_tpt = 1
                ,affected_period = 0.50
                ,n_clusters_affected = 1
                ,case='Case 4'
                ,sheetname = '1 Contaminated Small Balanced'
                )
execute_program(T = 10
                ,N = 27
                ,with_epidemic = TRUE
                ,change_param = c(0.10,0.20,0.30,0.40)
                ,start_tpt = 1
                ,affected_period = 0.50
                ,n_clusters =2
                ,n_clusters_affected = 1
                ,case='Case 4'
                ,sheetname = '1 Contaminated Small Unbalanced'
                )
execute_program(T = 10
                ,N = 31
                ,with_epidemic = TRUE
                ,change_param = c(0.10,0.20,0.30,0.40)
                ,start_tpt = 1
                ,affected_period = 0.50
                ,n_clusters=5
                ,n_clusters_affected = 1
                ,case='Case 4'
                ,sheetname = '1 Contaminated Small Unbalanced'
                )
execute_program(T = 50
                ,N = 50
                ,with_epidemic = TRUE
                ,change_param = c(0.10,0.20,0.30,0.40)
                ,start_tpt = 1
                ,affected_period = 0.50
                ,n_clusters_affected = 1
                ,case='Case 4'
                ,sheetname = '1 Contaminated Large Balanced'
                )


###all contaminated---------------
#### Case 1: contamination in all clusters, short period, no change in parameters ------
execute_program(T = 20
                ,N = 20
                ,with_epidemic = TRUE
                ,change_param = 0
                ,start_tpt = 1
                ,affected_period = 0.25
                ,n_clusters_affected = 'all'
                ,case='Case 1'
                ,sheetname = 'All Contaminated Small Balanced'
)
execute_program(T = 50
                ,N = 50
                ,with_epidemic = TRUE
                ,change_param = 0
                ,start_tpt = 1
                ,affected_period = 0.25
                ,n_clusters_affected = 'all'
                ,case='Case 1'
                ,sheetname = 'All Contaminated Large Balanced'
)
#### Case 2: contamination in all clusters, short period, with change in parameters ------
execute_program(T = 20
                ,N = 20
                ,with_epidemic = TRUE
                ,change_param = c(0.10,0.20,0.30,0.40)
                ,start_tpt = 1
                ,affected_period = 0.25
                ,n_clusters_affected = 'all'
                ,case='Case 2'
                ,sheetname = 'All Contaminated Small Balanced'
)
execute_program(T = 50
                ,N = 50
                ,with_epidemic = TRUE
                ,change_param = c(0.10,0.20,0.30,0.40)
                ,start_tpt = 1
                ,affected_period = 0.25
                ,n_clusters_affected = 'all'
                ,case='Case 2'
                ,sheetname = 'All Contaminated Large Balanced'
)
#### Case 3: contamination in all clusters, long period, no change in parameters ------
execute_program(T = 20
                ,N = 20
                ,with_epidemic = TRUE
                ,change_param = 0
                ,start_tpt = 1
                ,affected_period = 0.50
                ,n_clusters_affected = 'all'
                ,case='Case 3'
                ,sheetname = 'All Contaminated Small Balanced'
)
execute_program(T = 50
                ,N = 50
                ,with_epidemic = TRUE
                ,change_param = 0
                ,start_tpt = 1
                ,affected_period = 0.50
                ,n_clusters_affected = 'all'
                ,case='Case 3'
                ,sheetname = 'All Contaminated Large Balanced'
)
#### Case 4: contamination in all clusters, long period, with change in parameters ------
execute_program(T = 20
                ,N = 20
                ,with_epidemic = TRUE
                ,change_param = c(0.10,0.20,0.30,0.40)
                ,start_tpt = 1
                ,affected_period = 0.50
                ,n_clusters_affected = 'all'
                ,case='Case 4'
                ,sheetname = 'All Contaminated Small Balanced'
)
execute_program(T = 50
                ,N = 50
                ,with_epidemic = TRUE
                ,change_param = c(0.10,0.20,0.30,0.40)
                ,start_tpt = 1
                ,affected_period = 0.50
                ,n_clusters_affected = 'all'
                ,case='Case 4'
                ,sheetname = 'All Contaminated Large Balanced'
)



