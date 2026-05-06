### code Hendriek Boshuizen 
### Hendriek.Boshuizen@wur.nl
### Here we simulate data from a quadratic relationship between the exposure and the outcome (SBP), 
### simulate different measurement errors and
### apply different statistical methods to analyse these data

### the overall code for the paper used in the paper is divided in two parts:
### 1. making data for simulation and running the simulation (this code) plus storing the results on disc
### 2. visualizing results : data stored on file are read and used for visualisation

### This first part of the code comes in 3 variants: using a linear relation(this code), using a quadratic relation, and using a linear relation and then 
### dichotomize the outcome and apply logistic regression in stead of linear regression


### In this version Berkson error is calculation from error prone x using a prediction model
### fixed categories made small as otherwise lot of data
### fall outside outer categories
   
  require(foreach)
  library(doParallel)
  cl<-makeCluster(3, outfile = "debugquad.txt")
  registerDoParallel(cl)  
  
  rm(sim.dt)
  xtypes <-c("lognormal","normal","uniform")

  foreach(xx = 1:3) %dopar% {
##  for (xx in 1:1) {  ## in case no parallel computing is used
    require(data.table)
    require(reshape2)
    require(survival)
    require(gamlss)
   
    xtype <- xtypes[xx]
 
    
    ##############################################
    
    ##   simulate  datasets
    
    ##############################################
    
    
    
    #############################################################
    
    ##### simulation categorical variables and confounders ######
    
    #############################################################
    
    
    ## mean and std approximately as in duplo study
    nsim <- 500 # number of simulations per paralel process was 500
    nrec <- 500  # number of records in a single dataset
    distmean <- 100 # mean value of the exposure
    sigma <-
      0.24 # standard deviation of the exposure in the population on the logscale
    n <- nsim * nrec ## total number of simulated records
    
    # sd is approx 25 on the original scale
    
    
    ## simulate true exposure from a lognormal distribution
    if (xtype == "lognormal"){
      mu <- log(distmean) - 0.5 * sigma^2
      exptrue <- exp(rnorm(n, mu, sigma))
    }
    
    ## simulate true exposure from a normal distribution
    
    if (xtype == "normal"){
      exptrue <- rnorm(n, distmean, distmean*sqrt(exp(sigma^2)-1))
    }
    summary(exptrue)
    
    ## simulate true exposure from a uniform distribution
    
    if (xtype == "uniform"){
       exptrue <- distmean+distmean*sqrt(exp(sigma^2)-1)*(runif(n)-0.5)*sqrt(12) 
    }
    

    ## simulate another predictor (confounder1 and confounder2 are confounders, confounder0 not a confounder
    ## make them on average 0 with sd 1
    ## confounder0 is not a confounder because it is unrelated to intake but it is a determinant of the outcome
    
    sdexp <- distmean*sqrt(exp(sigma^2)-1)

    confounder0 <-      rnorm(n, 0, 1) 
    
    
    ## confounder1 that is positively associated with both outcome and exposure
    confounder1 <-
      sqrt(1-(0.02*sdexp)^2)*rnorm(n, 0, 1) + exptrue * 0.02 - mean(exptrue)* 0.02 ## confounder1 is confounder with same distribution as confounder0 and the same predictive value
    ## confounder2 of equal influence but correlated negatively with exposure
    confounder2 <-
      sqrt(1-(0.02*sdexp)^2)*rnorm(n, 0, 1) - exptrue * 0.02 + mean(exptrue)* 0.02 ## confounder2 is confounder with same distribution as confounder0 and the same predictive value
    
    ### changed 0.02 above to 0.033 to increasing the correlation from slightly over 0.5 to 0.9 approximately
    ###
    ### this was changed back to 0.02 because otherwise the residual method does not work very well
    ### as most of the variance of the exposure is taken out 
    
    ## print some summary information
    summary(cbind(confounder0, confounder1, confounder2)) 
    ## when exptrue is lognormal, confounder1 is not completely normal too
    cor(cbind(confounder0,confounder1,confounder2, exptrue))
  
    #############################################################
    
    ## Simulate the outcomes
    
    #############################################################
    
    ### simulate outcome SBP (true average systolic bloodpressure) in realistic range
    ### simulate sbp in range around 120-140
    ### intercept is 134 
    
    ### we want the range of SBP to be equal over all simulation (variance = 91 + 9 =100)
    ### however, exptrue and confounder0 are independent, while exptrue and confounder1/2
    ### are dependent. 
    ### this means that the variance of a*(exptrue-90)^2 + b*confounder will be different
    ### namely var(exptrue)*(2*a*(exptrue-90))^2 + var(conf)*b^2 +2 cov (2*a*b*(exptrue-90)) 
    ### we use the empirical values (variances and covariance calculated from the data) instead
    
    ### we use 5 as coefficient for the confounder, as the sd of (exp-90)^2  is approx
    ### 38 times as large as that of exp, here the coefficient of exp^2 is 25 times as small
    ### (so it is a slightly stronger effect). for a equally strong confounding effect
    ### we would need a coefficient of 4.56, but we made this 5
    ### this means variance from confounder is 25 and we need to add only 75
    
    ### SBP1: weak effect of exposure
    
    vartoadd <-  75-0.001^2*var((exptrue-90)^2) - 2*cov(5*confounder1,0.001*(exptrue-90)^2)
     
    SBP1<-132-  0.001* mean((exptrue-90)^2)+0.001*(exptrue-90)^2+5*(confounder0) +rnorm(n,0,sqrt(vartoadd))
    summary(SBP1)
    
    ### SBP2: strong effect of exposure
    
    vartoadd <-  75-0.004^2*var((exptrue-90)^2) - 2*cov(5*confounder1,0.004*(exptrue-90)^2)
    SBP2<-132-  0.004* mean((exptrue-90)^2)+0.004*(exptrue-90)^2+5*(confounder0) +rnorm(n,0,sqrt(vartoadd))
    summary(SBP2)
   
    ### SBP0: no effect of exposure
    
    vartoadd <-  75
    SBP0<-132+5*(confounder0)+rnorm(n,0,sqrt( vartoadd))
    summary(SBP0)
    
    ### similar outcomes with confounder 
    
    vartoadd <-  75-0.001^2*var((exptrue-90)^2) - 2*cov(5*confounder1,0.001*(exptrue-90)^2)
    SBP1C <- 132 -  0.001* mean((exptrue-90)^2)+0.001*(exptrue-90)^2+5*(confounder1) +rnorm(n,0,sqrt(vartoadd))
    summary(SBP1C)
    
    ### Strong effect with pos confounder
    vartoadd <-  75-0.004^2*var((exptrue-90)^2) - 2*cov(5*confounder1,0.004*(exptrue-90)^2)
    SBP2C <- 132 -  0.004* mean((exptrue-90)^2)+0.004*(exptrue-90)^2+5*(confounder1)+rnorm(n,0,sqrt(vartoadd))
    summary(SBP2C)
   
    ### null model
    vartoadd <-  75
    SBP0C <- 132+5*(confounder1)+rnorm(n,0,sqrt(vartoadd))
    summary(SBP0C)
    
    ### with a negative confounder 
    vartoadd <-  75-0.001^2*var((exptrue-90)^2) - 2*cov(5*confounder1,0.001*(exptrue-90)^2)
    SBP1C2<-132-  0.001* mean((exptrue-90)^2)+0.001*(exptrue-90)^2+5*(confounder2)+rnorm(n,0,sqrt( vartoadd))
    summary(SBP1C2)
     
    ### with a negative confounder and strong effect
    vartoadd <-  75-0.004^2*var((exptrue-90)^2) - 2*cov(5*confounder1,0.004*(exptrue-90)^2)
    SBP2C2<-132-  0.004* mean((exptrue-90)^2)+0.004*(exptrue-90)^2+5*(confounder2)+rnorm(n,0,sqrt(vartoadd))
    summary(SBP2C)
    
    ### null model
    
    ## with a confounder in a null model
    vartoadd <-  75
    SBP0C2 <-
      132 + 5 * (confounder2) + rnorm(n, 0, sqrt( vartoadd)) ## no effect of exposure and confounding
    summary(SBP0C2)
    
    ## print some summaries
    summary(cbind( SBP0, SBP1, SBP2))
    summary(cbind( SBP0C, SBP1C, SBP2C))
    summary(cbind( SBP0C2, SBP1C2, SBP2C2))
    
    cor(SBP2C2,confounder2 ) # 0.4 negative confounding gives same range of SBP 
    cor(SBP2C,confounder1 )  ### 0.68 and crude correlation is from confounder MINUS from exposure
    library(ppcor)
    ### print results in table 
    ### table in paper made with nsim 5000
    pcor(as.matrix(cbind(SBP2,confounder0,(exptrue-90)^2),ncol=3))$estimate
    pcor(as.matrix(cbind(SBP2C,confounder1,(exptrue-90)^2),ncol=3))$estimate
    pcor(as.matrix(cbind(SBP2C2,confounder2,(exptrue-90)^2),ncol=3))$estimate
    cor(cbind((exptrue-90)^2,confounder1 ,confounder2, confounder0))
    cor(cbind(exptrue,confounder1 ,confounder2, confounder0))
     
   
    #### simulate exposure with classical additive error
    
    x1 <- exptrue + rnorm(n, 0, 20) ## small error
    x2 <- exptrue + rnorm(n, 0, 50) ## large error
    x3 <-       pmax(x2, 0) ## truncated large error (x always positive) ## not used in paper
    
    ## print some correlations
    cor(x1, exptrue)
    cor(x2, exptrue)
    cor(exptrue,confounder1)
    cor(SBP2C,confounder1) # 2 is strong effect
    cor(SBP2C2,confounder2)
    cor(SBP2C,confounder2)
    cor(SBP2C2,confounder2)
    cor(SBP2C,exptrue)
    
    ## simulate a reference measurement to make calibrated exposure with
    
    reference <- exptrue + rnorm(n, 0, 20) ## small error
    
    ### simulate multiplicative error
    ### this will give a larger mean value for the x1m/x2m than exptrue
    ### However that will not influence results
    
    x1m <- exptrue * (1 + rnorm(n, 0, 0.2)) ## small error
    x2m <- exptrue * (1 + rnorm(n, 0, 0.5)) ## large error
    x3m <- pmax(x2m, 0)
    cor(x1m, exptrue)
    cor(x2m, exptrue)
    
    ### simulate x variables with Berkson error
    ### choice: small error : variable with Berson error 0.92 correlation
    ### large error 0.686 correlation
    
    x1b<- predict(lm(exptrue~pspline(x1,df=6)))  # gives exactly same correlation as between x1 and exptrue
    x2b<- predict(lm(exptrue~pspline(x2,df=6))) 
    x3b <-
      pmax(x2b, 0) # not needed as this will not easily be negative
    
    ### for uniform distribution at the edges residuals are not quit what they should be, as it is hard to realize
    ### berkson error at the ends of the distribution
    ### print some correlations
    cor(x1b,exptrue)
    cor(x1,exptrue)
    cor(x2b,exptrue)
    cor(x2,exptrue)
    cor(x2m,exptrue)
    summary(cbind(x1b,x1,x2b,x2,exptrue))
   
    cbind(var(x1b), var(exptrue))
    cbind(var(x2b), var(exptrue))
    cor(x1b, exptrue)
    cor(x2b, exptrue)
    
    
    ##########  collect the simulated data in a data table
    
    
    sim.dt <-
      data.table(
        SBP0,
        SBP1,
        SBP2,
        SBP0C,
        SBP1C,
        SBP2C,
        SBP0C2,
        SBP1C2,
        SBP2C2,
        x1,
        x2,
        x3,
        x1m,
        x2m,
        x3m,
        x1b,
        x2b,
        x3b,
        exptrue,
        reference,
        confounder0,
        confounder1,
        confounder2
      )
   
     ### names of exposure variables that will be used
    ### the variables with 3 in the name put negative values
    ### to zero. They are not used
    exposureNames <- c("exptrue", "x1", "x2",
                       "x1m", "x2m", "x1b", "x2b")
    
   
    ####################################################################################
    
    ####                    analyse the simulated datasets
    
    ####################################################################################
    
    ## we have now the following variants of data
    ## 3 association strengths exposure -outcome (incl. null)  ---- 3
    ## with / without confounder                               ---- 2
    ## 7 exposure (1 true, 6 with error)                       ---- 7
    ## quintiles and fixed categories                          ---- 2
    ##                                          total: 2x3x2x7 =  284
    
    
    ## we are using 10 methods of analysis:
    ##     - 1. naive: make categories and calculate difference with lowest category
    ##     - 2. calibrated: as naive, but using calibrated exposure (only differs when confounded)
    ## and results that can be compared to the true effect But different for fixed categories
    ##     - 3. use naive effect plus average exposure per category to calculate effect per unit of exposure
    ##     - 4. as above, but use average calibrated exposure
    ##     - 5. as above but use average reference exposure in uncalibrated categories
    ##     - 6. lastly, also a continuous exposure
    ##     - 7. continuous also calibrated
    ##     - 8. Methode ruth: divide by sqrt(att factor)
    ##     - 9. ditto but att factor without confounder adjustment
    ##     - 10. as 2 but with residual method
    ##     - 11. as 10 but per unit of exposure (stands to 9 as 3 to 2)
    ##     - 12  as 2 but using splines for calibration
    ##     - 13  as 12 but per unit of exposure
    ##     - 14 combining 10 and 12
    ##     - 15 combining 11 and 13
    ## for true exposure (error=0) we also do these, although they might not make sense\
    ## however, results can indicate the extra uncertainty from calibration
    
    
    ## create objects to store the output of the simulations (estimate beta per unit of exposure)
    ##### output: datatable with each row a result for the 4 quintiles
    ##### variables indicating what this is, namely: simno, confounder (y/n), strength association, exposure error, analysis type
    #### number of records = nsim*3 (conf) *3 (strength of association)*2 (cattypes)* 7 (#error) * 7 (methods) 
    #### effects of continuous analyses are stored in Q2
    nError <- 7  ## can be changed later
    nAnalyses <- 2 * 3 * 3 * (15 * nError) ## per simulated dataset
    NAf <- as.numeric(rep(NA, nsim * nAnalyses))
    results <-
      data.table(
        simno = NAf,
        conf = NAf,
        strength.association = NAf,
        error = NAf,
        analysis = NAf,
        categorytype = NAf,
        Q1 = NAf,
        ## Q1 will contain intercept
        Q2 = NAf,
        Q3 = NAf,
        Q4 = NAf,
        Q5 = NAf,
        cov2 = NAf,
        cov3 = NAf,
        cov4 = NAf,
        cov5 = NAf,
        T1 = NAf,
        T2 = NAf,
        ## true values, (no longer differences!)
        T3 = NAf,
        T4 = NAf,
        T5 = NAf,
        X1 = NAf,## mean value of category
        X2 = NAf,
        X3 = NAf,
        X4 = NAf,
        X5 = NAf,
        X1r = NAf,## mean value of category from reference measurement
        X2r = NAf,
        X3r = NAf,
        X4r = NAf,
        X5r = NAf
      )
    
    
    
    i<-1
    for (i in 1:nsim) {
      if (i %% 20 == 0)
        print(paste(xtype,i,Sys.time()))
      ### make a small data.table for use in this loop
      index <-
        ((i - 1) * nrec + 1):(i * nrec) ## pointer to records for current dataset
      loop.dt <- sim.dt[index,] ## current data
      
      ## repeat for options;
      confounder=0
      association=2
      cattype=1
      loopno = 0 ## number of the loop needed to store results at the right place
      for (confounder in c(0, 1, 2)) {
        for (association in c(0, 1, 2)) {
          for (cattype in 1:2) {
            loopno = loopno + 1
            ## fill data table with descriptive information
            ## in a single loop we we fill nError * 7 records
            index2 <-
              ((i - 1) * nAnalyses + (loopno - 1) * nError * 15 + 1):((i - 1) * nAnalyses +
                                                                       (loopno) * nError * 15) ## index in the result files
            results[index2, simno := i]
            results[index2, conf := confounder]
            results[index2, strength.association := association]
            results[index2, categorytype := cattype]
           
             ## get data for analysis from dataset
            
            SBP <- paste0("SBP", association)
            if (confounder == 1)
              SBP <-
              paste0(SBP, "C")           else  if (confounder == 2)
              SBP <-
              paste0(SBP, "C2")  ## SBP is current outcome in this loop
            
            if (confounder == 0)
              C <-
              "confounder0" else if (confounder == 1) 
              C <- "confounder1"    else
              C <- "confounder2" ## C is current confounder in this loop
            true.effect <- c(0, 0.025, 0.1)[association + 1]
            
            ### loop over every exposure with different error
            ### different errors are represented by different variable names in the datafile
            err <- 1
            for (err in 0:6) {
               
              X <- exposureNames[err + 1]
              ### make calibrated exposure (2= with splines)
              mod <-
                lm(paste0("reference~", X, " + ", C), data = loop.dt) ## make calibrated exposure
             
              data4gamlss <- data.frame(x=loop.dt[,get(X)], ref=loop.dt[,reference],c=loop.dt[,get(C)])
              mod2 <- gamlss(ref ~ c + pbm(x,mono="up"), data = data4gamlss,control = gamlss.control(trace=FALSE))
              cal <- predict(mod) # calibrated exposure
              calspline <- fitted(mod2) # calibrated exposure using splines
              attenuationfactor <- coef(mod)[2]
              mod.without.confounder <-
                lm(paste0("reference~", X ), data = loop.dt) 
              attenuationfactor2 <- coef(mod.without.confounder)[2]
              
              ## save calibrated exposure
              
              loop.dt[, cali := cal]
              loop.dt[, calispline := calspline]
             
             ### make calibrated exposure with residual method
              
              ### first make residuals
              makeres <-  lm(paste0(X, " ~ ", C), data = loop.dt) 
              resids <- residuals(makeres) + mean(loop.dt[,get(X)])  ### mean X is added to make the magnitude more natural
              ### for calibration we need a reference measurement for the residuals
              makeresref <-
                lm(paste0("reference ~ " , C), data = loop.dt) ## make calibrated exposure
              residsref <- residuals(makeresref) + mean(loop.dt[,reference])  ### mean X is added to make the magnitude more natural
              ### store them in the datatable
              loop.dt[, res := resids]
              loop.dt[, resref := residsref]
              cor(resids,residsref)
              
              ### as  res/resref are independent of C, adjusting for C might not be needed
              ### but it can also not do a lot of harm here
              ### calibration models. again 2 = using splines
              modres <-
                lm(paste0("resref ~ res + ", C), data = loop.dt) ## make calibrated exposure residual method
              cal10 <-predict(modres) 
              ### save in datatable
              loop.dt[, calresi := cal10] 
              
              data4gamlss$res <- loop.dt$res
              data4gamlss$resref <- loop.dt$resref
              modres2 <- gamlss(resref ~ pbm(res, mono="up") + c, data=data4gamlss,control = gamlss.control(trace=FALSE))
            
              cal14 <-fitted(modres2)
              loop.dt[, calresispline := cal14] 
              
              ## make quintiles both of original and calibrated variables
              loop.dt[, catx := cut(get(X),
                                    quantile(get(X), probs = seq(0, 1, .2)),
                                    include.lowest = TRUE)]
              loop.dt[, catcal := cut(cal, quantile(cal, probs = seq(0, 1, .2)), include.lowest = TRUE)]
              loop.dt[, catcalres := cut(cal10, quantile(cal10, probs = seq(0, 1, .2)), include.lowest = TRUE)]
              if ( length(unique(cal14)) >1 & abs(max(cal14)- min(cal14))>0.0001) loop.dt[, catcalresspline := cut(cal14, quantile(cal14, probs = seq(0, 1, .2)), include.lowest = TRUE)] else
                loop.dt[, catcalresspline := 1]
              loop.dt[, catcalspline := cut(calspline, quantile(calspline, probs = seq(0, 1, .2)), include.lowest = TRUE)]
              if (cattype == 2) {
                loop.dt[, catx := cut(get(X),
                                      c(-1000, 80, 95, 105, 120, 2000),
                                      include.lowest = TRUE)]
                loop.dt[, catcal := cut(cal,
                                        c(-1000, 80, 95, 105, 120, 2000),
                                        include.lowest = TRUE)]
                loop.dt[, catcalres := cut(cal10,
                                        c(-1000, 80, 95, 105, 120, 2000),
                                        include.lowest = TRUE)]
                loop.dt[, catcalspline := cut(calspline,
                                                 c(-1000, 80, 95, 105, 120, 2000),
                                                 include.lowest = TRUE)]
                if( length(unique(cal14)) >1  & abs(max(cal14)- min(cal14))>0.0001) loop.dt[, catcalresspline := cut(cal14,
                                           c(-1000, 80, 95, 105, 120, 2000),
                                           include.lowest = TRUE)] else
                                             loop.dt[, catcalresspline := 1]
                
              }
              
              
              summary(loop.dt)
              
              ## add the means of the calibrated measurement to categories of X and of calibrated X
              loop.dt[, meancat := mean(get(X)), by = catx] ## uncalibrated
              loop.dt[, meancatcal := mean(cali), by = catcal] ## calibrated means of calibrated categories
              loop.dt[, meancatcal2 := mean(reference), by = catx] ## reference means of uncalibrated categories
              loop.dt[, meancatcalres := mean(calresi), by = catcalres] ## reference means of uncalibrated categories
              loop.dt[, meancatcalresspline := mean(calresispline), by = catcalresspline] ## reference means of uncalibrated categories
              loop.dt[, meancatcalspline := mean(calispline), by = catcalspline] ## reference means of uncalibrated categories
              loop.dt[,catcalresspline]
              loop.dt[,meancatcalresspline]
  ## method 1: naive 
              meanx <- unique(loop.dt$meancat)
              meanx <- meanx[order(meanx)]
              
              if ((cattype == 2) & (min(loop.dt[,get(X)])>80)) meanx <- c(0,meanx)
              if ((cattype == 2) & (max(loop.dt[,get(X)])<=120)) meanx <- c(meanx,0)
   ## method 2: calibrated
              meanxcal <- unique(loop.dt$meancatcal)
              meanxcal <- meanxcal[order(meanxcal)]
  ## if there are no cases in the lowest/highest categorie this goes wrong           
              if ((cattype == 2) & (min(loop.dt$cali)>80)) meanxcal <- c(NA,meanxcal)
              if ((cattype == 2) & (max(loop.dt$cali)<=120)) meanxcal <- c(meanxcal,NA)
  ### method 5
              meanxcal2 <- unique(loop.dt$meancatcal2)
              meanxcal2 <- meanxcal2[order(meanxcal2)]
              if ((cattype == 2) & (min(loop.dt$reference)>80)) meanxcal2 <- c(NA,meanxcal2)
              if ((cattype == 2) & (max(loop.dt$reference)<=120)) meanxcal2 <- c(meanxcal2,NA)
  ### method 10            
              meanxcalres <- unique(loop.dt$meancatcalres)
              meanxcalres <- meanxcalres[order(meanxcalres)]
              ## if there are no cases in the lowest/highest categorie this goes wrong           
              if ((cattype == 2) & (min(loop.dt$calresi)>80)) meanxcalres <- c(NA,meanxcalres)
              if ((cattype == 2) & (min(loop.dt$calresi)>95)) meanxcalres <- c(NA,meanxcalres)
              if ((cattype == 2) & (max(loop.dt$calresi)<=120)) meanxcalres <- c(meanxcalres,NA)
              if ((cattype == 2) & (max(loop.dt$calresi)<=105)) meanxcalres <- c(meanxcalres,NA)
  
  ### method  12          
              meanxcalspline <- unique(loop.dt$meancatcalspline)
              meanxcalspline <- meanxcalspline[order(meanxcalspline)]
              ## if there are no cases in the lowest/highest categorie this goes wrong           
              if ((cattype == 2) & (min(loop.dt$calispline)>80)) meanxcalspline <- c(NA,meanxcalspline)
              if ((cattype == 2) & (min(loop.dt$calispline)>95)) meanxcalspline <- c(NA,meanxcalspline)
              if ((cattype == 2) & (max(loop.dt$calispline)<=120)) meanxcalspline <- c(meanxcalspline,NA)
              if ((cattype == 2) & (max(loop.dt$calispline)<=105)) meanxcalspline <- c(meanxcalspline,NA)
              
  ### method  14          
              meanxcalresspline <- unique(loop.dt$meancatcalresspline)
              meanxcalresspline <- meanxcalresspline[order(meanxcalresspline)]
              ## if there are no cases in the lowest/highest categorie this goes wrong           
              if ((cattype == 2) & (min(loop.dt$calresispline)>80)) meanxcalresspline <- c(NA,meanxcalresspline)
              if ((cattype == 2) & (min(loop.dt$calresispline)>95)) meanxcalresspline <- c(NA,meanxcalresspline)
              if ((cattype == 2) & (max(loop.dt$calresispline)<=120)) meanxcalresspline <- c(meanxcalresspline,NA)
              if ((cattype == 2) & (max(loop.dt$calresispline)<=105)) meanxcalresspline <- c(meanxcalresspline,NA)
              
  ### method 5 separate references measurements
              meanxcal2 <- unique(loop.dt$meancatcal2)
              meanxcal2 <- meanx[order(meanxcal2)]
              if ((cattype == 2) & (min(loop.dt[,get(X)])>80)) meanxcal2 <- c(NA,meanxcal2)
              if ((cattype == 2) & (max(loop.dt[,get(X)])<=120)) meanxcal2 <- c(meanxcal2,NA)
             
                
  ### run the models (naive and on calibrated data)
              modobs.error <-
                summary(lm(paste(SBP, "~ catx +", C), data = loop.dt))[[4]] ## coefficients
              modobs.cal <-
                summary(lm(paste(SBP, "~ catcal +", C), data = loop.dt))[[4]] ## coefficients
              if ( sum(!is.na(meanxcalres))>1)  modobs.calres <-
                summary(lm(paste(SBP, "~ catcalres + ", C), data = loop.dt))[[4]] 
              if ( sum(!is.na(meanxcalspline))>1)  modobs.calspline <-
                summary(lm(paste(SBP, "~ catcalspline + ", C), data = loop.dt))[[4]] ## coefficients
              if ( sum(!is.na(meanxcalresspline))>1)  modobs.calresspline <-
                summary(lm(paste(SBP, "~ catcalresspline + ", C), data = loop.dt))[[4]] ## coefficients
              
              ## store the coefficients
              ## start with adding descriptives to results file
              indexStart <-
                index2[err * 15 + 1] ## index to store first result
              results[indexStart:(indexStart + 14), error := err]
              results[indexStart:(indexStart + 14), analysis := 1:15]
              
      ### naive        
              results[indexStart, c("Q1", "Q2", "Q3", "Q4", "Q5") := as.list(modobs.error[1:5, 1])]
              
              results[indexStart, c("X1", "X2", "X3", "X4", "X5") := as.list(meanx)]
              
               
      ### method 8
              results[indexStart + 7, c( "Q2", "Q3", "Q4", "Q5") := as.list(modobs.error[2:5, 1]/sqrt(attenuationfactor))]
      ### method 9
              results[indexStart + 8, c( "Q2", "Q3", "Q4", "Q5") := as.list(modobs.error[2:5, 1]/sqrt(attenuationfactor2))]
      ### method 2
              results[indexStart+1, c("X1", "X2", "X3", "X4", "X5") := as.list(meanxcal)]
              
              results[indexStart+1, c("X1r", "X2r", "X3r", "X4r", "X5r") := as.list(meanxcal2)]
      ### method 10
              results[indexStart+9, c("X1", "X2", "X3", "X4", "X5") := as.list(meanxcalres)]
       
      ### method 12
              results[indexStart + 11, c("X1", "X2", "X3", "X4", "X5") := as.list(meanxcalspline)]
      ### method 14
              if (sum(!is.na(meanxcalresspline))>1) results[indexStart + 13, c("X1", "X2", "X3", "X4", "X5") := as.list(meanxcalresspline)]
              
                     
      ### method 2        
              ### for fixed categories this only works if data are available, 
              ### here we only make data in the case all are present, or the lower or upper category is missing
              if (length(modobs.cal[,1])==6) results[indexStart + 1, c("Q1", "Q2", "Q3", "Q4", "Q5") := as.list(modobs.cal[1:5, 1])]
              if (length(modobs.cal[,1])==5 & is.na(meanxcal[1]))  results[indexStart + 1, c("Q2", "Q3", "Q4", "Q5") := as.list(modobs.cal[1:4, 1])]
              if (length(modobs.cal[,1])==5 & is.na(meanxcal[5]))  results[indexStart + 1, c("Q1","Q2", "Q3", "Q4") := as.list(modobs.cal[1:4, 1])]
              if (length(modobs.cal[,1])==4 & is.na(meanxcal[1]) & is.na(meanxcal[5]))  results[indexStart + 1, c("Q2", "Q3", "Q4") := as.list(modobs.cal[1:3, 1])]
      
        ### method 10       
              ### for fixed categories this only works if data are available, 
              ### here we only make data in the case all are present, or the lower or upper category is missing
              if (length(modobs.calres[,1])==6) results[indexStart + 9, c("Q1", "Q2", "Q3", "Q4", "Q5") := as.list(modobs.calres[1:5, 1])]
              if (length(modobs.calres[,1])==5 & is.na(meanxcalres[1]))  results[indexStart + 9, c("Q2", "Q3", "Q4", "Q5") := as.list(modobs.calres[1:4, 1])]
              if (length(modobs.calres[,1])==5 & is.na(meanxcalres[5]))  results[indexStart + 9, c("Q1","Q2", "Q3", "Q4") := as.list(modobs.calres[1:4, 1])]
              if (length(modobs.calres[,1])==4 & is.na(meanxcalres[1]) & is.na(meanxcalres[5]))  results[indexStart + 9, c("Q2", "Q3", "Q4") := as.list(modobs.calres[1:3, 1])]
              
        ### NB to calculate the true effect we need to find the mean difference
        ### of true exposure and multiply this with the effect.
              
        ### method 12 
              
              ### for fixed categories this only works if data are available, 
              ### here we only make data in the case all are present, or the lower or upper category is missing
              if (length(modobs.calspline[,1])==6) results[indexStart + 11, c("Q1", "Q2", "Q3", "Q4", "Q5") := as.list(modobs.calspline[1:5, 1])]
              if (length(modobs.calspline[,1])==5 & is.na(meanxcalspline[1]))  results[indexStart + 11, c("Q2", "Q3", "Q4", "Q5") := as.list(modobs.calspline[1:4, 1])]
              if (length(modobs.calspline[,1])==5 & is.na(meanxcalspline[5]))  results[indexStart + 11, c("Q1","Q2", "Q3", "Q4") := as.list(modobs.calspline[1:4, 1])]
              if (length(modobs.calspline[,1])==4 & is.na(meanxcalspline[1]) & is.na(meanxcalspline[5]))  results[indexStart + 11, c("Q2", "Q3", "Q4") := as.list(modobs.calspline[1:3, 1])]
              
        ### method 14       
              ### for fixed categories this only works if data are available, 
              ### here we only make data in the case all are present, or the lower or upper category is missing
              if (length(modobs.calresspline[,1])==6) results[indexStart + 13, c("Q1", "Q2", "Q3", "Q4", "Q5") := as.list(modobs.calresspline[1:5, 1])]
              if (length(modobs.calresspline[,1])==5 & is.na(meanxcalresspline[1]))  results[indexStart + 13, c("Q2", "Q3", "Q4", "Q5") := as.list(modobs.calresspline[1:4, 1])]
              if (length(modobs.calresspline[,1])==5 & is.na(meanxcalresspline[5]))  results[indexStart + 13, c("Q1","Q2", "Q3", "Q4") := as.list(modobs.calresspline[1:4, 1])]
              if (length(modobs.calresspline[,1])==4 & is.na(meanxcalresspline[1]) & is.na(meanxcalresspline[5]))  results[indexStart + 13, c("Q2", "Q3", "Q4") := as.list(modobs.calresspline[1:3, 1])]
              
              
              loop.dt[, cattrue := cut(exptrue,
                                       quantile(exptrue, probs = seq(0, 1, .2)),
                                       include.lowest = TRUE)]
              if (cattype == 2)
                loop.dt[, cattrue := cut(exptrue,
                                         c(-1000, 80, 95, 105, 120, 2000),
                                         include.lowest = TRUE)]
    ### this has no place in the file, is also redundant because same as error = 0 
              meantrue <-
                unlist(unique(loop.dt[, mean(exptrue), by = cattrue])[order(cattrue)][, 2])
              if ((cattype == 2) & (min(loop.dt[,exptrue])>80)) meantrue <- c(NA,meantrue)
              if ((cattype == 2) & (max(loop.dt[,exptrue])<=120)) meantrue <- c(meantrue,NA)
              meancat1 <- unlist(meantrue[1])
              
              
              ##### change here for linear /quadratic !!!!
              
             true.effect.adj <- c(132, 120, 112)[association + 1] -
                c(0, 0.025, 0.1)[association + 1] * (meantrue[1:5])
              
              
              ### save this in results file to be used for method 1 and 2 (+10,12,14)
              results[indexStart, c("T1", "T2", "T3", "T4", "T5") := as.list(c(
                true.effect.adj[1],
                true.effect.adj[2],
                true.effect.adj[3],
                true.effect.adj[4],
                true.effect.adj[5]
              ))]
              results[indexStart + 1, c("T1", "T2", "T3", "T4", "T5") :=as.list( c(
                true.effect.adj[1],
                true.effect.adj[2],
                true.effect.adj[3],
                true.effect.adj[4],
                true.effect.adj[5]
              ))]
              
               results[indexStart + 9, c("T1", "T2", "T3", "T4", "T5") :=as.list( c(
                true.effect.adj[1],
                true.effect.adj[2],
                true.effect.adj[3],
                true.effect.adj[4],
                true.effect.adj[5]
              ))]
               
               results[indexStart + 11, c("T1", "T2", "T3", "T4", "T5") :=as.list( c(
                 true.effect.adj[1],
                 true.effect.adj[2],
                 true.effect.adj[3],
                 true.effect.adj[4],
                 true.effect.adj[5]
               ))]
               
               results[indexStart + 13, c("T1", "T2", "T3", "T4", "T5") :=as.list( c(
                 true.effect.adj[1],
                 true.effect.adj[2],
                 true.effect.adj[3],
                 true.effect.adj[4],
                 true.effect.adj[5]
               ))]
              
              ### make differences with first category
              ### does not work when no cases in lowest category
              true.effect.adj <- true.effect.adj - meancat1
              true.effect.adj <- true.effect.adj[-1]
              ### method 1
              lowlim1 <-
                (modobs.error[2:5, 1] - 1.96 * modobs.error[2:5, 2])
              highlim1 <-
                (modobs.error[2:5, 1] + 1.96 * modobs.error[2:5, 2])
          ###  method 2-->3
              if (length(modobs.cal[,1])==6) lowlim2 <-
                (modobs.cal[2:5, 1] - 1.96 * modobs.cal[2:5, 2])
              if (length(modobs.cal[,1])==6) highlim2 <-
                (modobs.cal[2:5, 1] + 1.96 * modobs.cal[2:5, 2])
          ### method 11 (from 10)
              if (length(modobs.calres[,1])==6) lowlim10 <-
                (modobs.calres[2:5, 1] - 1.96 * modobs.calres[2:5, 2])
              if (length(modobs.calres[,1])==6) highlim10 <-
                (modobs.calres[2:5, 1] + 1.96 * modobs.calres[2:5, 2])
              
          ### method 13 (from 12)
              if (length(modobs.calspline[,1])==6) lowlim12 <-
                (modobs.calspline[2:5, 1] - 1.96 * modobs.calspline[2:5, 2])
              if (length(modobs.calspline[,1])==6) highlim12 <-
                (modobs.calspline[2:5, 1] + 1.96 * modobs.calspline[2:5, 2])
              
          ### method 15 (from 14)
              if (length(modobs.calresspline[,1])==6) lowlim14 <-
                (modobs.calresspline[2:5, 1] - 1.96 * modobs.calresspline[2:5, 2])
              if (length(modobs.calresspline[,1])==6) highlim14 <-
                (modobs.calresspline[2:5, 1] + 1.96 * modobs.calresspline[2:5, 2])
              
       ### coverage of true effect       
              results[indexStart, c("cov2", "cov3", "cov4", "cov5") := as.list(lowlim1 <
                                                                                 true.effect.adj &
                                                                                 highlim1 > true.effect.adj)]
              results[indexStart + 1, c("cov2", "cov3", "cov4", "cov5") := as.list(lowlim2 <
                                                                                     true.effect.adj &
                                                                                     highlim2 > true.effect.adj)]
                                                                                      
              
              
              ### analysis of the SBP increase per unit of exposure as derived from the categorical analysis
              meanvalues.error <-
                unique(loop.dt[, meancat, by = catx])[order(catx)][1:5, meancat]
              meanvalues.cal <-
                unique(loop.dt[, meancatcal, by = catcal])[order(catcal)][1:5, meancatcal]
              meanvalues.cal2 <-
                unique(loop.dt[, meancatcal2, by = catx])[order(catx)][1:5, meancatcal2]
          ### method 11
              meanvalues.calres <-
              unique(loop.dt[, meancatcalres, by = catcalres])[order(catcalres)][1:5, meancatcalres]
          ### method 13
              meanvalues.calspline <-
                unique(loop.dt[, meancatcalspline, by = catcalspline])[order(catcalspline)][1:5, meancatcalspline]
          ### method 15
              meanvalues.calresspline <-
                unique(loop.dt[, meancatcalresspline, by = catcalresspline])[order(catcalresspline)][1:5, meancatcalresspline]
              
               
              delta.error <-
                meanvalues.error[2:5] - meanvalues.error[1]
              delta.cal <- meanvalues.cal[2:5] - meanvalues.cal[1]
              delta.cal2 <- meanvalues.cal2[2:5] - meanvalues.cal2[1]
              delta.calres <- meanvalues.calres[2:5] - meanvalues.calres[1]
              delta.calspline <- meanvalues.calspline[2:5] - meanvalues.calspline[1]
              delta.calresspline <- meanvalues.calresspline[2:5] - meanvalues.calresspline[1]
              
              dif.error <- modobs.error[2:5] / delta.error
              if (length(modobs.cal[,1])==6) dif.cal <- modobs.cal[2:5] / delta.cal
              dif.cal2 <- modobs.error[2:5] / delta.cal2
              if (length(modobs.calres[,1])==6) dif.calres <- modobs.calres[2:5] / delta.calres
              if (length(modobs.calspline[,1])==6) dif.calspline <- modobs.calspline[2:5] / delta.calspline
              if (length(modobs.calresspline[,1])==6) dif.calresspline <- modobs.calresspline[2:5] / delta.calresspline
              
               
              results[indexStart + 2, c("Q2", "Q3", "Q4", "Q5") := as.list(dif.error)]
              results[indexStart + 3, c("Q2", "Q3", "Q4", "Q5") := as.list(dif.cal)]
              results[indexStart + 4, c("Q2", "Q3", "Q4", "Q5") := as.list(dif.cal2)]
              results[indexStart + 10, c("Q2", "Q3", "Q4", "Q5") := as.list(dif.calres)]
              results[indexStart + 12, c("Q2", "Q3", "Q4", "Q5") := as.list(dif.calspline)]
              results[indexStart + 14, c("Q2", "Q3", "Q4", "Q5") := as.list(dif.calresspline)]
              
              
              ## coverage (excluding effect of calibration): true values are 0, 0.025 and 0.1
              lowlim.error <-
                (modobs.error[2:5, 1] - 1.96 * modobs.error[2:5, 2]) / delta.error
              highlim.error <-
                (modobs.error[2:5, 1] + 1.96 * modobs.error[2:5, 2]) / delta.error
              if (length(modobs.cal[,1])==6) lowlim.cal <-
                (modobs.cal[2:5, 1] - 1.96 * modobs.cal[2:5, 2]) / delta.cal
              if (length(modobs.cal[,1])==6) highlim.cal <-
                (modobs.cal[2:5, 1] + 1.96 * modobs.cal[2:5, 2]) / delta.cal
              lowlim.cal2 <-
                (modobs.error[2:5, 1] - 1.96 * modobs.error[2:5, 2]) / delta.cal2
              highlim.cal2 <-
                (modobs.error[2:5, 1] + 1.96 * modobs.error[2:5, 2]) / delta.cal2
         ### method 11
              if (length(modobs.calres[,1])==6) lowlim.calres <-
                (modobs.calres[2:5, 1] - 1.96 * modobs.calres[2:5, 2]) / delta.calres
              if (length(modobs.calres[,1])==6) highlim.calres <-
                (modobs.calres[2:5, 1] + 1.96 * modobs.calres[2:5, 2]) / delta.calres
         ### method 13
              if (length(modobs.calspline[,1])==6) lowlim.calspline <-
                (modobs.calspline[2:5, 1] - 1.96 * modobs.calspline[2:5, 2]) / delta.calspline
              if (length(modobs.calspline[,1])==6) highlim.calspline <-
                (modobs.calspline[2:5, 1] + 1.96 * modobs.calspline[2:5, 2]) / delta.calspline
          ### method 15
              if (length(modobs.calresspline[,1])==6) lowlim.calresspline <-
                (modobs.calresspline[2:5, 1] - 1.96 * modobs.calresspline[2:5, 2]) / delta.calresspline
              if (length(modobs.calresspline[,1])==6) highlim.calresspline <-
                (modobs.calresspline[2:5, 1] + 1.96 * modobs.calresspline[2:5, 2]) / delta.calresspline
              
              
              
              results[indexStart + 2, c("cov2", "cov3", "cov4", "cov5") := as.list(lowlim.error <
                                                                                     true.effect &
                                                                                     highlim.error > true.effect)]
              results[indexStart + 3, c("cov2", "cov3", "cov4", "cov5") := as.list(lowlim.cal <
                                                                                     true.effect &
                                                                                     highlim.cal > true.effect)]
              results[indexStart + 4, c("cov2", "cov3", "cov4", "cov5") := as.list(lowlim.cal2 <
                                                                                     true.effect &
                                                                                     highlim.cal2 > true.effect)]
              
              results[indexStart + 10, c("cov2", "cov3", "cov4", "cov5") := as.list(lowlim.calres <
                                                                                     true.effect &
                                                                                     highlim.calres > true.effect)]
              results[indexStart + 12, c("cov2", "cov3", "cov4", "cov5") := as.list(lowlim.calspline <
                                                                                      true.effect &
                                                                                      highlim.calspline > true.effect)]
              results[indexStart + 14, c("cov2", "cov3", "cov4", "cov5") := as.list(lowlim.calresspline <
                                                                                      true.effect &
                                                                                      highlim.calresspline > true.effect)]
              
              
              ## continous : true values are 0, 0.025 and 0.1
              
              mod.cont.error <-
                summary(lm(paste(SBP, "~", X, " + ", C), data = loop.dt))[[4]]
              mod.cont.cal <-
                summary(lm(paste(SBP, "~ cali + ", C), data = loop.dt))[[4]]
              
              
              results[indexStart + 5, Q2 := mod.cont.error[2, 1]]
              results[indexStart + 6, Q2 := mod.cont.cal[2, 1]]
              
              lowlim3 <-
                (mod.cont.error[2, 1] - 1.96 * mod.cont.error[2, 2])
              highlim3 <-
                (mod.cont.error[2, 1] + 1.96 * mod.cont.error[2, 2])
              
              results[indexStart + 5, cov2 := (lowlim3 < true.effect &
                                                 highlim3 > true.effect)]
              
              lowlim4 <-
                (mod.cont.cal[2, 1] - 1.96 * mod.cont.cal[2, 2])
              highlim4 <-
                (mod.cont.cal[2, 1] + 1.96 * mod.cont.cal[2, 2])
              
              results[indexStart + 6, cov2 := (lowlim4 < true.effect &
                                                 highlim4 > true.effect)]
              
              ## also store the true effects to calculate rmse later
              true.effect <- unlist(true.effect)
              results[indexStart + 2, c("T1", "T2", "T3", "T4", "T5") :=
                        list(true.effect[1],
                             true.effect[1],
                             true.effect[1],
                             true.effect[1],
                             true.effect[1])]
              results[indexStart + 3, c("T1", "T2", "T3", "T4", "T5") :=
                        list(true.effect[1],
                             true.effect[1],
                             true.effect[1],
                             true.effect[1],
                             true.effect[1])]
              results[indexStart + 4, c("T1", "T2", "T3", "T4", "T5") :=
                        list(true.effect[1],
                             true.effect[1],
                             true.effect[1],
                             true.effect[1],
                             true.effect[1])]
              results[indexStart + 5, c("T1", "T2", "T3", "T4", "T5") :=
                        list(true.effect[1],
                             true.effect[1],
                             true.effect[1],
                             true.effect[1],
                             true.effect[1])]
              results[indexStart + 6, c("T1", "T2", "T3", "T4", "T5") :=
                        list(true.effect[1],
                             true.effect[1],
                             true.effect[1],
                             true.effect[1],
                             true.effect[1])]
              
              results[indexStart + 10, c("T1", "T2", "T3", "T4", "T5") :=
                        list(true.effect[1],
                             true.effect[1],
                             true.effect[1],
                             true.effect[1],
                             true.effect[1])]
              results[indexStart + 12, c("T1", "T2", "T3", "T4", "T5") :=
                        list(true.effect[1],
                             true.effect[1],
                             true.effect[1],
                             true.effect[1],
                             true.effect[1])]
              
              
              results[indexStart + 14, c("T1", "T2", "T3", "T4", "T5") :=
                        list(true.effect[1],
                             true.effect[1],
                             true.effect[1],
                             true.effect[1],
                             true.effect[1])]
              
              
              
              
              
            }
            
          }
        }
        
      }
    }
    
    warnings()
    
    #check results
    results[, mean(Q2,  na.rm = TRUE), by = list(conf, strength.association, error, analysis)]
    print(summary(results))
    
    ### save in different objects depending on xtype
    
    save(file = paste0("resultssim_quad9", xtype,".RData"), results)
    
   
    
  }
  
  stopImplicitCluster()
  
