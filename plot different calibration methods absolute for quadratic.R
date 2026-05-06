################################################

## plotting the results of the simulations
## for quintiles categories and the quadratic model

## parameter main=true makes separate plots for the main part of the paper
## main = false makes a pdf as supplementary material
##

## the program reads three input files (made by the code "simulation and analysis quadratic model.R") 
## names of these files are : resultssim_quad9normal.RData, resultssim_quad9lognormal.RData, resultssim_quad9uniform.RData
## in this example they are stored in the subdirectory "Rscripts and data last versions"

## plots produced are stored in the subdirectory "plots"
## change directory names when needed

################################################

require(ggplot2)
require(grid)
require(gridExtra)
require(data.table)

main=FALSE

###############################
###
### plotting parameters
###
###############################

errorstrength <- 1

###############################################################################
#
# read in data to be plotted and assign identifying labels to print on the plot
#
###############################################################################

resultlist <- list(3)
load("Rscripts and data last versions/resultssim_quad9normal.RData")
resultlist[[1]] <- results
load("Rscripts and data last versions/resultssim_quad9lognormal.RData")
resultlist[[2]] <- results
load("Rscripts and data last versions/resultssim_quad9uniform.RData")
resultlist[[3]] <- results
simlabel <- "quadratic"
labels.distribution<- c( "protein:normal dist.","protein: lognormal dist.","protein: uniform dist.")

if (!main) pdf("plots/Supplement_S3_quad.pdf")

#############################################################################

## prepare data for plots 

#############################################################################


summary(results$Q2)
dist=1

p <- vector("list", 9) # list to store the plots
#for (dist in 1:3){
for (dist in 1:3){
  
### make a dataframe containing the true relation between intake and SBP
 
  sigma=0.24
  distmean=100
  if (dist==2) { # (xtype == "lognormal")
    mu <- log(distmean) - 0.5 * sigma^2
    intake <- exp(qnorm(seq(0.005,0.995,0.01), mu, sigma))
  }
  
  if (dist == 1 ) { ## "normal"
    
    intake <- qnorm(seq(0.005,0.995,0.01),distmean, distmean*(exp(sigma)-1))
  }
  
  ### SBP2: strong effect of exposure
  ### SBP is simulated with:
  ### SBP2C <- 132 -  0.004* mean((exptrue-90)^2) + 0.004*(exptrue-90)^2 + 5*(confounder1) + rnorm(n,0,sqrt(vartoadd))
  ### We make a data.frame containing the simulated relation 
  
  ### start with making a list with intake values
  
  if (dist==3) { #} (xtype == "uniform")
    
    intake <- distmean+distmean*(exp(sigma)-1)*(seq(-0.495,0.495,0.01))*sqrt(12) 
  }
  
  ### calculate the true SBP for this intake values
 
  SBPtrue <- 132-  0.004* mean((intake-90)^2)+0.004*(intake-90)^2
  
  ###  make a dataframe called linedat to store the results 
  ###  store this for each combination of error-type and analysis-type
  linedat <- data.frame(SBP=rep(SBPtrue,15), intake=rep(intake,15),  
                       errortype =rep( c("additive", "multiplicative", "Berkson"), each=5*length(intake)),
                      analysis = rep(rep(c(1,2,10,12,14),each=length(intake)),3) )                  

  ### calculate the outcome to be plotted
  
  ### make absolute results by adding the SBP of the reference category (Q1)
  results<- resultlist[[dist]]
  results$absQ1<-results$Q1
  results$absQ2<-results$Q2+results$Q1
  results$absQ3<-results$Q3+results$Q1
  results$absQ4<-results$Q4+results$Q1
  results$absQ5<-results$Q5+results$Q1
  ### make difference in exposure (D) by subtracting the exposure in the reference category
  results$D2<-results$X2-results$X1
  results$D3<-results$X3-results$X1
  results$D4<-results$X4-results$X1
  results$D5<-results$X5-results$X1

  ## print some checks 
  summary(subset(results,analysis <3 |  analysis==12 )) 
  nrow(subset(results,analysis <3|  analysis==12 )) # 126000 /500*2 
  
  ### 126 different datapoint per simulation
  ### 7 errortypes * 3 assoc * 2 cattypes * 3 confounding
  
  ### make a dataframe in the format needed to plot using ggplot
  
  plotdata <-
    melt(results,
         id = c("simno", "conf", "strength.association","categorytype",
                "error", "analysis"))

  ## add variables to characterize each data point
  plotdata$type <- "NA"
  plotdata[grep("Q", plotdata$variable), ]$type <- "relative effect"
  plotdata[grep("absQ", plotdata$variable), ]$type <- "absolute effect"
  plotdata[grep("Q1", plotdata$variable), ]$type <- "intercept"
  plotdata[grep("absQ1", plotdata$variable), ]$type <- "absolute effect"
  plotdata[grep("T", plotdata$variable), ]$type <- "true effect"
  plotdata[grep("cov", plotdata$variable), ]$type <- "coverage"
  plotdata[grep("X", plotdata$variable), ]$type <- "X"
  plotdata[grep("X[0-9]r", plotdata$variable), ]$type <- "Xr"
  plotdata[grep("D", plotdata$variable), ]$type <- "dX"
  table(plotdata$type, plotdata$variable)  
  plotdata$errortype<-NA
  plotdata$errortype[(plotdata$error==0 )] <- "none"
  plotdata$errortype[(plotdata$error %in% 1:2 )] <- "additive"
  plotdata$errortype[(plotdata$error %in% 3:4 )] <- "multiplicative"
  plotdata$errortype[(plotdata$error %in% 5:6 )] <- "Berkson"
  plotdata$error2 <- 2-plotdata$error%%2 
 
  ### put the continuous data in a separate category
  ### do in two steps because somehow changing the variable messes up the comparison (was OK in earlier R version)
  ### select only the data that are existing (where continues data exist)
  isQcont<- plotdata$variable == "Q2" & plotdata$analysis > 5 & plotdata$analysis <8
  isCovcont <- plotdata$variable == "cov2" & plotdata$analysis > 5 & plotdata$analysis <8
  
  plotdata$variable <- factor(plotdata$variable, levels=c(levels(plotdata$variable), "contin."))
  plotdata[ isQcont|isCovcont, ]$variable <- "contin."
  
  ## give the continuous data same analysis type as comparable categorical analyses
  plotdata[plotdata$analysis == 6, ]$analysis <- 3 
  plotdata[plotdata$analysis == 7, ]$analysis <- 4
  
  ## clean up the dataset by removing redundant (NA) values
  plotdata <- plotdata[!is.na(plotdata$value), ] ## remove missings

  ## select only X data that we are going to plot and put those in pl1
  pl1<-subset(plotdata,type=="X" & (analysis<3 | analysis==10| analysis==12 | analysis==14)) # only naive and calibrated 
  ## renaming of variables 
  pl1$X<-pl1$value
  pl1$variable<-substr(pl1$variable,2,2)
  ## select only Y data that we are going to plot and put those in pl0
  pl0<-subset(plotdata,(type=="absolute effect") & (analysis<3 | analysis==10| analysis==12 | analysis==14))
  pl0$variable<-substr(pl0$variable,5,5)
  pl0$SBP<-pl0$value
  ## merge them
  plotdatamerged <- merge(subset(pl0,select = -c(value,type)),subset(pl1,select = -c(value,type)))

  ### In Q5 there are missing values 318 times dus to lacking data in those categories. this is relatively few so we ignore 
  ### The reference values not used here, but can be added
  ### by also making and merging pl2 and then making Xr equal to X for method == 1
  
 
  ## make labels for plots
  analysis.labs = c("naive", "calibrated", "cal. residuals","spline cal.","spline cal. resid.")
  names(analysis.labs) <- c("1", "2", "10","12","14")
  error.labs = c("no error", "mod. error (additive)", "strong error (additief)","mod. error (multipl)", "strong error (multipl)",
                 "mod. error (Berkson)", "strong error (Berkson)")
  errortype.labs = c("additive","multiplicative","Berksonian")
  names(error.labs) <- c("0", "1", "2","3","4","5","6")
  error2.labs<- c("1" = "moderate error" , "2" = "large error")
  lab.confounding <- c("no confounder", "weak confounder", "strong confounder")
  reg.coef <- c(0,0.025,0.1)

  
  ##################################################################
  
  #####  loop for plotting
  
  ##################################################################
  
  for (association in 2:2 ) for (confounder in 0:2)    {
 
    lab.assoc<-c(" no association"," moderate association", " strong association" )[1+association]
    lab.confounding<-c(" no confounding"," pos. confounding"," neg. confounding") [confounder+1]
    
    ### plot absolute effects for the 5 quintiles 
    
    ### make plotdf which contains results of simulations of variables with error
     
    plotdf<-subset(
      plotdatamerged, error2 ==errorstrength  &
        strength.association == association &
        conf == confounder  & error > 0   & categorytype == 1)

    ### plotdf2 contains results of simulations of error free variables
    plotdf2<-subset(
      plotdatamerged, 
      strength.association == association &
        conf == confounder & error == 0 & categorytype == 1
    )
    
    ## plotdf2 needs to be in every plot so it is copied 3 times
    plotdf2.1 <-plotdf2
    plotdf2.2 <-plotdf2
    plotdf2.3 <-plotdf2
    
    plotdf2.1$errortype <- "additive" 
    plotdf2.2$errortype <- "multiplicative"
    plotdf2.3$errortype <- "Berkson"
    plotdf2<-data.table( rbind(plotdf2.1,plotdf2.2,plotdf2.3))
    plotdf <- data.table(plotdf)
 
    ## calculate the mean values of X and Y per category over the 500 simulated datasets for the errorfree data
    Q1value <- plotdf2[variable == "1", .(Q1value = mean(SBP)), by = .(errortype,analysis)] # for variable with error
    X1value <- plotdf2[variable == "1", .(X1value = mean(X)), by = .(errortype,analysis)] # for variable with error
    Q2value <- plotdf2[variable == "2", .(Q2value = mean(SBP)), by = .(errortype,analysis)] # for variable with error
    X2value <- plotdf2[variable == "2", .(X2value = mean(X)), by = .(errortype,analysis)] # for variable with error
    Q3value <- plotdf2[variable == "3", .(Q3value = mean(SBP)), by = .(errortype,analysis)] # for variable with error
    X3value <- plotdf2[variable == "3", .(X3value = mean(X)), by = .(errortype,analysis)] # for variable with error
    Q4value <- plotdf2[variable == "4", .(Q4value = mean(SBP)), by = .(errortype,analysis)] # for variable with error
    X4value <- plotdf2[variable == "4", .(X4value = mean(X)), by = .(errortype,analysis)] # for variable with error
    Q5value <- plotdf2[variable == "5", .(Q5value = mean(SBP)), by = .(errortype,analysis)] # for variable with error
    X5value <- plotdf2[variable == "5", .(X5value = mean(X)), by = .(errortype,analysis)] # for variable with error
    
    ## and for the continuous data
    Qvalues <- plotdf[, .(count = .N, SBP = mean (SBP), X=mean(X) , sdSBP = sd(SBP), sdX=sd(X)), by=.(variable, errortype, analysis) ]
    Qvalues2 <- plotdf2[, .(count = .N, SBP = mean (SBP), X=mean(X) , sdSBP = sd(SBP), sdX=sd(X)), by=.(variable, errortype, analysis) ]
    Qvalues$seSBP = with(Qvalues, sdSBP/sqrt(count))
    Qvalues2$seSBP = with(Qvalues2, sdSBP/sqrt(count))
    Qvalues$seX = with(Qvalues, sdX/sqrt(count))
    Qvalues2$seX = with(Qvalues2, sdX/sqrt(count))
    
    linedat1Q <- merge(Q1value,X1value,by= c("errortype","analysis"))
    linedat2Q <- merge(Q2value,X2value,by= c("errortype","analysis"))
    linedat3Q <- merge(Q3value,X3value,by= c("errortype","analysis"))
    linedat4Q <- merge(Q4value,X4value,by= c("errortype","analysis"))
    linedat5Q <- merge(Q5value,X5value,by= c("errortype","analysis"))
    names(linedat1Q) <-  names(linedat2Q) <- names(linedat3Q) <-
      names(linedat4Q) <- names(linedat5Q) <-c("errortype", "analysis" ,"SBP","intake")
    linedatQ1 <- rbind(linedat1Q,linedat2Q,linedat3Q,linedat4Q,linedat5Q)
    
    ## repeat for the variables with error
    
    plotdf <- data.table(plotdf)
    Q1value2 <- plotdf[variable == "1", .(Q1value = mean(SBP)), by = .(errortype,analysis)] # for variable with error
    X1value2 <- plotdf[variable == "1", .(X1value = mean(X)), by = .(errortype,analysis)] # for variable with error
    Q2value2 <- plotdf[variable == "2", .(Q2value = mean(SBP)), by = .(errortype,analysis)] # for variable with error
    X2value2 <- plotdf[variable == "2", .(X2value = mean(X)), by = .(errortype,analysis)] # for variable with error
    Q3value2 <- plotdf[variable == "3", .(Q3value = mean(SBP)), by = .(errortype,analysis)] # for variable with error
    X3value2 <- plotdf[variable == "3", .(X1value = mean(X)), by = .(errortype,analysis)] # for variable with error
    Q4value2 <- plotdf[variable == "4", .(Q4value = mean(SBP)), by = .(errortype,analysis)] # for variable with error
    X4value2 <- plotdf[variable == "4", .(X4value = mean(X)), by = .(errortype,analysis)] # for variable with error
    Q5value2 <- plotdf[variable == "5", .(Q5value = mean(SBP)), by = .(errortype,analysis)] # for variable with error
    X5value2 <- plotdf[variable == "5", .(X5value = mean(X)), by = .(errortype,analysis)] # for variable with error
    
    linedat1Q2 <- merge(Q1value2,X1value2,by= c("errortype","analysis"))
    linedat2Q2 <- merge(Q2value2,X2value2,by= c("errortype","analysis"))
    linedat3Q2 <- merge(Q3value2,X3value2,by= c("errortype","analysis"))
    linedat4Q2 <- merge(Q4value2,X4value2,by= c("errortype","analysis"))
    linedat5Q2 <- merge(Q5value2,X5value2,by= c("errortype","analysis"))
    names(linedat1Q2) <-  names(linedat2Q2) <- names(linedat3Q2) <-
    names(linedat4Q2) <- names(linedat5Q2) <-c("errortype", "analysis" ,"SBP","intake")
    
    ## combinae in overall dataframe
    linedatQ2 <- rbind(linedat1Q2,linedat2Q2,linedat3Q2,linedat4Q2,linedat5Q2)
    linedat <-subset(linedat, SBP<145) ## only plot for SBP values <145 
    
    
    
    ###################################################################
    #
    #  Select less or more plots
    #
    ###################################################################
    
    
    ###### to select only a few plots use the code below (also adjust facet in plot) 
    # plotdf <-subset(plotdf,errortype=="multiplicative" & analysis %in% c(1,10))
    # plotdf2 <-subset(plotdf2,errortype=="multiplicative" & analysis %in% c(1,10))
    # linedat <-subset(linedat,errortype=="multiplicative" & analysis %in% c(1,10))
    if (main){
      plotdf <-subset(plotdf, analysis %in% c(1,10) & errortype != "Berkson")
      plotdf2 <-subset(plotdf2, analysis %in% c(1,10) & errortype != "Berkson")
      linedat <-subset(linedat, analysis %in% c(1,10) & errortype != "Berkson" & SBP<140)
      linedatQ1 <-subset(linedatQ1, analysis %in% c(1,10) & errortype != "Berkson" & SBP<140)
      linedatQ2 <-subset(linedatQ2, analysis %in% c(1,10) & errortype != "Berkson"  & SBP<140)
      Qvalues <- subset(Qvalues, analysis %in% c(1,10) & errortype != "Berkson")
      Qvalues2 <- subset(Qvalues2, analysis %in% c(1,10) & errortype != "Berkson")
      
   
    }
    # adapt the label size depending on the number of plots 
    # this was optimized for the paper, needs tweeking if different number of plots is selected above
    if (main) lsize <- 11 else lsize <-8
    if (!main) subtitletag <- "blue: simulated relation; black: no error; red = with error" else subtitletag <- ""
    
    p[[(dist-1)*3+confounder+1]] <-
      ggplot(
        Qvalues,
        aes(y = SBP, x = X, group=variable)
      ) +
      ggtitle(paste0(labels.distribution[dist] ,"; ",lab.confounding)) +
      geom_point(data=Qvalues2, colour= "black", size=3, shape=16)+ 
      geom_point(data=Qvalues, colour="red",size=3, shape=16)+
      geom_errorbar(data=Qvalues2,aes(ymin=SBP-seSBP, ymax=SBP+seSBP), colour= "black", width=.2)+ 
      geom_errorbar(data=Qvalues, aes(ymin=SBP-seSBP, ymax=SBP+seSBP),colour="red", width=.2)+
      geom_line(data=linedat,aes(x=intake, y=SBP, group=1), colour= "blue",size=1) +
      geom_line(data=linedatQ1,aes(x=intake, y=SBP, group=1), colour= "black",size=1) +
      geom_line(data=linedatQ2,aes(x=intake, y=SBP, group=1), colour= "red",size=1) +
      labs(y = "SBP (mmHg)", x = "protein intake (g/d)") +
      facet_grid(
        errortype~ analysis ,
        labeller =labeller(errortype = label_value, analysis = analysis.labs))+
                theme_bw()+
      theme(plot.title = element_text(size=13), 
            plot.subtitle=element_text(size=10),
              axis.title.x = element_text(size = 13),
              axis.text.x = element_text(size = lsize+2),
              axis.text.y = element_text(size = 13),
              axis.title.y = element_text(size = 13),
              strip.text.x = element_text(size = lsize+2),
              strip.text.y = element_text(size = 12))+
      theme(panel.spacing.y = unit(1, "lines")) +
      scale_y_continuous(breaks = seq(125, 140, by = 5), limits = c(125, 140)) +
      scale_x_continuous(breaks = seq(50, 150, by = 50), limits = c(40, 160))
    ##  code used when only plotting a single errortype 
    ##  facet_wrap(~ analysis,
    ##  labeller =labeller(analysis = analysis.labs))
    
    ## plot result to pdf for supplementary material
    if (!main) print(p[[(dist-1)*3+confounder+1]])
    
    
    
  }
  
}

## plot to png for main paper
if (main) {
  p1 <- grid.arrange(p[[1]],p[[4]],p[[2]],p[[5]], nrow = 2)
  ggsave(file=paste0("plots/figure 4.png"),plot=p1,device="png",width=21, height=18,unit="cm")
  
  
}

if (!main) dev.off()
