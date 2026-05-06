################################################
## plotting the results of the simulations
## for fixed categories compared to first group
## in this version with quintiles
##
## parameter main=true makes separate plots for the main part of the paper
## main = false makes a pdf as supplementary material
##
## the program reads three input files (made by the code "simulation and analysis linear model binary version.R") 
## names of these files are : resultssim_lin9normal.RData, resultssim_lin9lognormal.RData, resultssim_lin9uniform.RData
## in this example they are stored in the subdirectory "Rscripts and data last versions"
##
## plots and rtf tables produced are stored in the subdirectory "plots"
## change directory names when needed

################################################

require(ggplot2)
require(data.table)
library(rtf)
require(grid)
require(gridExtra)

###############################
###
### plotting parameters
###
###############################

errorstrength <- 2
main <- FALSE

###############################################################################
#
# read in data to be plotted and assign identifying labels to print on the plot
#
###############################################################################

resultlist <- list(3)
load("Rscripts and data last versions/resultssim_lin9normal.RData")
resultlist[[1]] <- results
load("Rscripts and data last versions/resultssim_lin9lognormal.RData")
resultlist[[2]] <- results
load("Rscripts and data last versions/resultssim_lin9uniform.RData")
resultlist[[3]] <- results
simlabel <- "linear"
if (!main) pdf("plots/Supplement_S1_linear.pdf")
if (!main) rtf<- RTF(file="plots/outputtab_linear_all.rtf", width=8.5,height=11,font.size=10,omi=c(1,1,1,1))

#############################################################################

##   plot results

#############################################################################
labels.distribution<- c( "protein:normal dist.","protein: lognormal dist.","protein: uniform dist.")

#############################################################################

## prepare data for plots 

#############################################################################

summary(results$Q2)
dist=3
p <- vector("list", 9) ## save the plots made
names(results)

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
 
  if (dist==3) { #} (xtype == "uniform")
    
    intake <- distmean+distmean*(exp(sigma)-1)*(seq(-0.495,0.495,0.01))*sqrt(12) 
  }
  
  
  ### SBP2: strong effect of exposure
  ### simulated relationship
  ###  SBP2C2 <-   134 - 0.1*mean(exptrue) + 0.1 * exptrue +3 * (confounder2) + rnorm(n, 0, sqrt(vartoadd) ## strong effect of exposure and confounding
  
  ### We plot these between 50 and 150 for confounder=0
  
  SBPtrue <- 134 -  0.1* mean(intake)+0.1*intake
 
  ###  make a dataframe called linedatT to store the results 
  ###  store this for each combination of error-type and analysis-type
  linedatT<-data.frame(SBP=rep(SBPtrue,15), intake=rep(intake,15),  
                       errortype =rep( c("additive", "multiplicative", "Berkson"), each=5*length(intake)),
                       analysis = rep(rep(c(1,2,10,12,14),each=length(intake)),3) )
  
  ### make absolute results by adding the SBP of the reference category (Q1)
  results<- resultlist[[dist]]
  results$absQ1<-results$Q1
  results$absQ2<-results$Q2+results$Q1
  results$absQ3<-results$Q3+results$Q1
  results$absQ4<-results$Q4+results$Q1
  results$absQ5<-results$Q5+results$Q1
  results$D2<-results$X2-results$X1
  results$D3<-results$X3-results$X1
  results$D4<-results$X4-results$X1
  results$D5<-results$X5-results$X1
  summary(subset(results,analysis <3 | analysis==10 | analysis==12 | analysis==14)) # Q1 60 NA for normal Q5 15
  nrow(subset(results,analysis <3| analysis==10 | analysis==12 | analysis==14)) # 126000 /500*2 
 
   # old code selecting different analyses
  # subset(resultlist[[1]],is.na(Q1) & (analysis<3 | analysis==10 | analysis==12 | analysis==14))$error
  
  ## melt to format needed in ggplot
  plotdata <-
    melt(results,
         id = c("simno", "conf", "strength.association","categorytype",
                "error", "analysis"))
  unique(plotdata$variable)
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
  with(plotdata,table(error2,errortype))
  
  ## put the continuous data in a separate category
  ### do in two steps because somehow changing the variable messes up the comparison (was OK in earlier R version)
  isQcont<- plotdata$variable == "Q2" & plotdata$analysis > 5 & plotdata$analysis <8
  isCovcont <- plotdata$variable == "cov2" & plotdata$analysis > 5 & plotdata$analysis <8
  
  plotdata$variable <- factor(plotdata$variable, levels=c(levels(plotdata$variable), "contin."))
  plotdata[isQcont, ]$variable <- "contin."
  plotdata[isCovcont, ]$variable <- "contin."
  ## give the continuous data same analysis type as comparable categorical analyses
  plotdata[plotdata$analysis == 6, ]$analysis <- 3 
  plotdata[plotdata$analysis == 7, ]$analysis <- 4
 
  ## clean up the dataset by removing redundant (NA) values
  plotdata <- plotdata[!is.na(plotdata$value), ] ## remove missings
  ## select only X data that we are going to plot and put those in pl1
  pl1<-subset(plotdata,type=="X" & (analysis<3 | analysis==10| analysis==12 | analysis==14)) # only naive and calibrated analyses
  pl1$X<-pl1$value
  pl1$variable<-substr(pl1$variable,2,2)
  
  ## select only Y data that we are going to plot and put those in pl0
  pl0<-subset(plotdata,(type=="absolute effect") & (analysis<3 | analysis==10| analysis==12 | analysis==14))
  unique(pl0$variable)
  pl0$variable<-substr(pl0$variable,5,5)
  pl0$SBP<-pl0$value
  ## merge them
  plotdatamerged <- merge(subset(pl0,select = -c(value,type)),subset(pl1,select = -c(value,type)))
  
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
  association=2
  
  confounder=2
  
  method=1
  # hulp is introduced to make it easy to plot only a subset of the analyses
  # this can be use to single out one or more plots  
  # example:   
  # hulp <-(subset(
  #    plotdatamerged, error2 ==1 & 
  #      strength.association == 1 &
  #      conf == 2  & error ==3   & categorytype == 1 & analysis==14))
  
  
  
  for (association in 2:2 ) for (confounder in 0:2)    {
    ##################################################################
    
    ###  start plot loop ###
    
    ##################################################################
    
    lab.assoc<-c(" no association"," moderate association", " strong association" )[1+association]
    lab.confounding<-c(" no confounding"," pos. confounding"," neg. confounding") [confounder+1]
    
    ### plot absolute effects for the 5 quintiles 
    
    ### plotdf contains results of simulations of variables with error
    plotdf<-subset(
      plotdatamerged, error2 ==errorstrength &
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
   
    Q1value <- plotdf2[variable == "1", .(Q1value = mean(SBP)), by = .(errortype,analysis)] # for variable with error
    X1value <- plotdf2[variable == "1", .(X1value = mean(X)), by = .(errortype,analysis)] # for variable with error
    Q5value <- plotdf2[variable == "5", .(Q5value = mean(SBP)), by = .(errortype,analysis)] # for variable with error
    X5value <- plotdf2[variable == "5", .(X5value = mean(X)), by = .(errortype,analysis)] # for variable with error

    Qvalues <- plotdf[, .(count = .N, SBP = mean (SBP), X=mean(X) , sdSBP = sd(SBP), sdX=sd(X)), by=.(variable, errortype, analysis) ]
    Qvalues2 <- plotdf2[, .(count = .N, SBP = mean (SBP), X=mean(X) , sdSBP = sd(SBP), sdX=sd(X)), by=.(variable, errortype, analysis) ]
    Qvalues$seSBP = with(Qvalues, sdSBP/sqrt(count))
    Qvalues2$seSBP = with(Qvalues2, sdSBP/sqrt(count))
    Qvalues$seX = with(Qvalues, sdX/sqrt(count))
    Qvalues2$seX = with(Qvalues2, sdX/sqrt(count))
    
    ## make a dataframe containing the errorfree 1st and 5th quintile results
    linedat1 <- merge(Q1value,X1value,by= c("errortype","analysis"))
    linedat5 <- merge(Q5value,X5value,by= c("errortype","analysis"))
    linedat <- merge(linedat1,linedat5,by= c("errortype","analysis"))
    linedat$slope <- with(linedat,(Q5value-Q1value) /(X5value-X1value))  # so we do it in the oldfashioned way
    linedat$intercept <- with(linedat,Q1value-slope*X1value)
    
    plotdf <- data.table(plotdf)
    Q1value2 <- plotdf[variable == "1", .(Q1value = mean(SBP)), by = .(errortype,analysis)] # for variable with error
    X1value2 <- plotdf[variable == "1", .(X1value = mean(X)), by = .(errortype,analysis)] # for variable with error
    Q5value2 <- plotdf[variable == "5", .(Q5value = mean(SBP)), by = .(errortype,analysis)] # for variable with error
    X5value2 <- plotdf[variable == "5", .(X5value = mean(X)), by = .(errortype,analysis)] # for variable with error
    
    linedat2.1 <- merge(Q1value2,X1value2,by= c("errortype","analysis"))
    linedat2.5 <- merge(Q5value2,X5value2,by= c("errortype","analysis"))
    linedat2 <- merge(linedat2.1,linedat2.5,by= c("errortype","analysis"))
    #linedat <- linedat[,slope := (Q4value - Q2value)/2,] # does not work for unknown reasons
    linedat2$slope <- with(linedat2,(Q5value-Q1value) /(X5value-X1value))  # so we do it in the oldfashioned way
    linedat2$intercept <- with(linedat2,Q1value-slope*X1value)
    
    ###################################################################################
    ##
    ##     to select only a few plots adapt the code below
    ##
    ###################################################################################
    
    
    if (main) {
    plotdf <-subset(plotdf, analysis %in% c(1,10) & errortype != "Berkson")
    plotdf2 <-subset(plotdf2, analysis %in% c(1,10) & errortype != "Berkson")
    linedat <-subset(linedat, analysis %in% c(1,10) & errortype != "Berkson")
    linedat2 <-subset(linedat2, analysis %in% c(1,10) & errortype != "Berkson")
    linedatT <-subset(linedatT, analysis %in% c(1,10) & errortype != "Berkson")
    Qvalues <- subset(Qvalues, analysis %in% c(1,10) & errortype != "Berkson")
    Qvalues2 <- subset(Qvalues2, analysis %in% c(1,10) & errortype != "Berkson")
    
   }
    # adapt the label size depending on the number of plots 
    # this was optimized for the paper, needs tweeking if different number of plots is selected above
    
    if (main) lsize <- 11 else lsize <-8
    
    p[[(dist-1)*3+confounder+1]] <-
      ggplot(data=Qvalues
        ,
        aes(y = SBP, x = X)
      ) +
      geom_abline(intercept=134-100*reg.coef[association+1], slope=reg.coef[association+1], colour= "blue", size=1) +
      geom_abline(data=linedat,aes(intercept=intercept, slope=slope), colour= "black", size=1) +
      geom_abline(data=linedat2,aes(intercept=intercept, slope=slope), colour= "red", size=1) +
      geom_line(data=linedatT,aes(x=intake, y=SBP, group=1), colour= "blue") +
      ggtitle(paste0(labels.distribution[dist] ,"; ",lab.confounding) ) +
      geom_point(data=Qvalues2, colour= "black", size=3, shape=16)+ 
      geom_point(data=Qvalues, colour="red",size=3, shape=16)+
      geom_errorbar(data=Qvalues2,aes(ymin=SBP-seSBP, ymax=SBP+seSBP), colour= "black", width=.2)+ 
      geom_errorbar(data=Qvalues, aes(ymin=SBP-seSBP, ymax=SBP+seSBP),colour="red", width=.2)+
      
      labs(y = "SBP (mmHg)", x = "protein intake (g/d)") +
      facet_grid(
        errortype ~ analysis ,
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
      theme(panel.spacing.y = unit(1, "lines"))+ 
      scale_y_continuous(breaks = seq(130, 140, by = 5),limits = c(130, 140))

    ## plot result to pdf for supplementary material
    if (!main) print(p[[(dist-1)*3+confounder+1]])
    ## print table with result to rtf 
    if(!main){  
      addHeader(rtf,title=paste(" distribution ",labels.distribution[dist]," confounder ",confounder))
      addText(rtf,"v3= slope between Q1 and Q5; v4= slope between Q2 and Q4")
      Q1value <- plotdf[variable == "1", .(Q1value = SBP), by = .(analysis,errortype)] # for variable with error
      X1value <- plotdf[variable == "1", .(X1value = X), by = .(analysis,errortype)] # for variable with error
      Q5value <- plotdf[variable == "5", .(Q5value = SBP), by = .(analysis,errortype)] # for variable with error
      X5value <- plotdf[variable == "5", .(X5value = X), by = .(analysis,errortype)] # for variable with error
      Q2value <- plotdf[variable == "2", .(Q2value = SBP), by = .(analysis,errortype)] # for variable with error
      X2value <- plotdf[variable == "2", .(X2value = X), by = .(analysis,errortype)] # for variable with error
      Q4value <- plotdf[variable == "4", .(Q4value = SBP), by = .(analysis,errortype)] # for variable with error
      X4value <- plotdf[variable == "4", .(X4value = X), by = .(analysis,errortype)] # for variable with error
      
      m1 <-merge(Q5value[, ID := .I],Q1value[, ID := .I])
      m2 <- merge(X5value[, ID := .I],X1value[, ID := .I])
      slope1 <- merge(m1,m2)
      slope1 <- slope1[,slope51 := (Q5value-Q1value)/(X5value-X1value)]
      m21 <- merge(Q4value[, ID := .I],Q2value[, ID := .I])
      m22 <- merge(X4value[, ID := .I],X2value[, ID := .I])
      slope2 <- merge(m21,m22)
      slope2 <- slope2[,slope42 := (Q4value-Q2value)/(X4value-X2value)]
      # for variable with error
      slop1 <-slope1[, .(mslope51 = mean(slope51)), by = .(analysis,errortype)] 
      slop2 <-slope2[, .(mslope42 = mean(slope42)), by = .(analysis,errortype)] 
      
      slope <- merge(slop1,slop2)
      setorder(slope, cols = "errortype")   
      addTable(rtf,as.data.frame(slope[,.(analysis,errortype,round(mslope51,4),round(mslope42,4))]))
      
      Q1value.0 <- plotdf2[variable == "1", .(Q1value = SBP), by = .(analysis,errortype)] # for variable without error
      X1value.0 <- plotdf2[variable == "1", .(X1value = X), by = .(analysis,errortype)] # for variable without error
      Q5value.0 <- plotdf2[variable == "5", .(Q5value = SBP), by = .(analysis,errortype)] # for variable without error
      X5value.0 <- plotdf2[variable == "5", .(X5value = X), by = .(analysis,errortype)] # for variable without error
      Q2value.0 <- plotdf2[variable == "2", .(Q2value = SBP), by = .(analysis,errortype)] # for variable without error
      X2value.0 <- plotdf2[variable == "2", .(X2value = X), by = .(analysis,errortype)] # for variable without error
      Q4value.0 <- plotdf2[variable == "4", .(Q4value = SBP), by = .(analysis,errortype)] # for variable without error
      X4value.0 <- plotdf2[variable == "4", .(X4value = X), by = .(analysis,errortype)] # for variable without error
      
      
      m1 <- merge(Q5value.0[, ID := .I],Q1value.0[, ID := .I])
      m2 <- merge(X5value.0[, ID := .I],X1value.0[, ID := .I])
      slope1.0 <- merge(m1,m2)
      slope1.0 <- slope1.0[,slope51 := (Q5value-Q1value)/(X5value-X1value)]
      m1 <- merge(Q4value.0[, ID := .I],Q2value.0[, ID := .I])
      m2 <- merge(X4value.0[, ID := .I],X2value.0[, ID := .I])
      slope2.0 <- merge(m1,m2)
      slope2.0 <- slope2.0[,slope42 := (Q4value-Q2value)/(X4value-X2value)]
      
      print( slope1.0[, .(mslope51 = mean(slope51)), by = .(analysis,errortype)] )# for variable with error
      print( slope2.0[, .(mslope42 = mean(slope42)), by = .(analysis,errortype)] ) 
      
    }
    
    
  }
  
}

## plot to png for main paper
if (main) {
  p1 <- grid.arrange(p[[1]],p[[4]],p[[2]],p[[5]], nrow = 2)
  ggsave(file=paste0("plots/figure 1.png"),plot=p1,device="png",width=21, height=18,unit="cm")
 }

if (!main) done(rtf)
if (!main) dev.off()

