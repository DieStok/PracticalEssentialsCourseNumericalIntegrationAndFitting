setwd("YOURWORKINGDIRECTORYHERE")
source("YOURGRINDSCRIPTHERE")
set.seed(12345)









##Solutions ordered by the new question format. (which runs only until 4A)









# The 3 experiments are indexed by their date:
exptA <- "2015-11-18"  #Experiment A
exptB <- "2015-12-14"  #Experiment B
exptC <- "2016-01-06"  #Experiment C



# First read and plot all data (for exptA):

fig3RA <- read.csv(paste("Fig3/",exptA,"_R.csv",sep="")) # Red
fig3GA <- read.csv(paste("Fig3/",exptA,"_G.csv",sep="")) # Green
plot(fig3RA$Time, fig3RA$OD, ylim=c(0,0.8),col="red",pch=".",xlab="Time (hr)",ylab="OD")
points(fig3GA$Time, fig3GA$OD, ylim=c(0,0.8),col="green",pch=".")

fig4 <- read.csv(paste("Fig4/",exptA,"_RG.csv",sep=""))
plot(fig4$Time, fig4$OD, ylim=c(0,0.8), col="blue", pch=".", xlab="Time (hr)", ylab="Total OD")

fig5 <- read.csv(paste("Fig5/flow_df_",exptA,".csv",sep=""))
fig5G <- fig5[fig5$Strain=="Green",]
fig5R <- fig5[fig5$Strain=="Red",]
plot(fig5G$time, fig5G$freq_mean, ylim=c(0,1), col="green",xlab="Time (hr)", ylab="Frequency")
points(fig5R$time, fig5R$freq_mean, col="red")

# Here the fitting starts 

model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    a <- q0/(q0+exp(-m*t))
    dN <- r*a*N*(1-(N/K)^v)
    return(list(dN))  
  }) 
}

#########################################
##Question 2A                          ##
#########################################


s <- c(N=0.124)
p <- c(K=0.6,r=0.4,m=2,q0=0.005,v=2)
free <- c("N",names(p))
lower <- rep(0, length(free)); lower[match("v",free)] <- 1; lower  # set lower bounds

#red
data3RA <- as.data.frame(cbind(fig3RA$Time,fig3RA$OD)); names(data3RA) <- c("time","N")
fit3RA <- fit(data3RA,free=free,fun=log,lower=lower,pch=".",legend=FALSE,tstep=0.1,main="Fig 3 a red monoculture (fun=log)",
              ymin = 0, ymax = 1)
summary(fit3RA)

#green
data3GA <- as.data.frame(cbind(fig3GA$Time,fig3GA$OD)); names(data3GA) <- c("time","N")
fit3GA <- fit(data3GA,free=free,fun=log,lower=lower,pch=".",legend=FALSE,tstep=0.1,main="Fig 3 a green monoculture (fun=log)",
              ymin = 0, ymax = 1)
summary(fit3GA)

#fitting without taking the log of the data:

fit3GANonLog <- fit(data3GA,free=free,lower=lower,pch=".",legend=FALSE,tstep=0.1,main="Fig 3 A green nonlog",
                    ymin = 0, ymax = 1)
summary(fit3GANonLog)
fit3RANonLog <- fit(data3RA,free=free,lower=lower,pch=".",legend=FALSE,tstep=0.1,main="Fig 3 A red nonlog",
                    ymin = 0, ymax = 1)
summary(fit3RANonLog)

hist(fit3GA$residuals, breaks = 20)
hist(fit3GANonLog$residuals, breaks = 20)


#Nouja, punt is: als je de log neemt dan zorg je dat discrepanties van kleine waardes (die je squared als je
#gaat fitten om de residuals te berekenen, en die daardoor kleiner worden) even zwaar worden meegeteld als van grotere waardes


########################################################################################
##                                                                                    ##
## Loop to read in data for figure 3A, B, and C, and store in list with associated fit##
##                                                                                    ##
########################################################################################

#I do this so I can easily calculate fits using Fig4A, and fig3A data later on, same for B, C.
expSeq = c(exptA, exptB, exptC)
names(expSeq) = c("exptA", "exptB", "exptC")
listFitsFigure3 <- list()

for (exp in seq_along(expSeq)) {
  
  dataFig3RLoop = read.csv(paste("Fig3/",expSeq[exp],"_R.csv",sep="")) # Red
  dataFig3GLoop = read.csv(paste("Fig3/",expSeq[exp],"_G.csv",sep="")) # Green
  
  sLoop <- c(N=0.124)
  pLoop <- c(K=0.6,r=0.4,m=2,q0=0.005,v=2)
  freeLoop <- c("N",names(pLoop))
  lowerLoop <- rep(0, length(freeLoop))
  lowerLoop[match("v",freeLoop)] <- 1 # set lower bounds
  
  #Red data and fit
  
  dataFig3RLoop <- as.data.frame(cbind(dataFig3RLoop$Time,dataFig3RLoop$OD))
  names(dataFig3RLoop) <- c("time","N")
  print(paste0("Fitting figure 3 data for the red bacteria for: ", names(expSeq[exp])))
  fitFig3RLoop <- fit(dataFig3RLoop,free=freeLoop,fun=log,lower=lowerLoop,
                   pch=".",legend=FALSE,tstep=0.1,main=paste0("Fig 3 Red for ",names(expSeq[exp])),
                   ymin = 0, ymax = 1, odes = model, parms = pLoop, state = sLoop)
  
  #Green data and fit
  
  dataFig3GLoop <- as.data.frame(cbind(dataFig3GLoop$Time,dataFig3GLoop$OD))
  names(dataFig3GLoop) <- c("time","N")
  print(paste0("Fitting figure 3 data for the green bacteria for: ", names(expSeq[exp])))
  fitFig3GLoop <- fit(dataFig3GLoop,free=freeLoop,fun=log,lower=lowerLoop,
                   pch=".",legend=FALSE,tstep=0.1,main=paste0("Fig 3 Green for ",names(expSeq[exp])),
                   ymin = 0, ymax = 1, odes = model, parms = pLoop, state = sLoop)
  
  #add to list under correct name
  listFitsFigure3[[paste0("Data", names(expSeq)[exp], "Fig3G")]] <- dataFig3GLoop
  listFitsFigure3[[paste0("Data", names(expSeq)[exp], "Fig3R")]] <- dataFig3RLoop
  listFitsFigure3[[paste0("Fit", names(expSeq)[exp], "Fig3G")]]  <- fitFig3GLoop
  listFitsFigure3[[paste0("Fit", names(expSeq)[exp], "Fig3R")]]  <- fitFig3RLoop
}




######################################################
# Question 2 b Taking the lower bound of v to be zero#
######################################################

lowerVZero <- rep(0, length(free)); lowerVZero[match("v",free)] <- 0; lowerVZero  # set lower bounds
data3RAVZero <- as.data.frame(cbind(fig3RA$Time,fig3RA$OD)); names(data3RAVZero) <- c("time","N")
fit3RAVZero <- fit(data3RAVZero,free=free,fun=log,lower=lowerVZero,pch=".",legend=FALSE,tstep=0.1,main="red fig3a v lower bound 0 (get very low v)",
                   ymin = 0, ymax = 1)
summary(fit3RAVZero)
#Completely unrealistic relation growth speed and pop size, see:
#normal:
curve(1-(x/0.6)^1, ylab = "Per capita growth if v = 1")
#ridiculous:
curve(1-(x/0.6)^0.26e-03, ylab = "Per capita growth if v is ridiculously low (0.26e-03), which the fit produces")
#in a very narrow range there is density-dependence
curve(1-(x/0.6)^0.26e-03, from = 0.124, ylab = "Per capita growth if v is ridiculously low, zoom in from starting OD of bacteria (on x)")

#Can fit figure 3 red, but cannot get correlation summary values, because V is not a square numeric matrix
#and the system is singular. Means we are fitting too many parameters: can't distinguish between v and r.
#predicted v = 0.26e-03, which is very small.

data3GAVZero <- as.data.frame(cbind(fig3GA$Time,fig3GA$OD)); names(data3GAVZero) <- c("time","N")
fit3GAVZero <- fit(data3GAVZero,free=free,fun=log,lower=lowerVZero,pch=".",legend=FALSE,tstep=0.1,main="green fig3a V lower bound 0",
                   ymin = 0, ymax = 1)
summary(fit3GAVZero)

#can fit figure 3 Green. Values barely change. 


###################################
#Exercise 2c fitting experiment 3B#
###################################

#What happens when lag phase not ignored? (Already the case here, since all parameters are free)


fig3RB <- read.csv(paste("Fig3/",exptB,"_R.csv",sep="")) # Red
fig3GB <- read.csv(paste("Fig3/",exptB,"_G.csv",sep="")) # Green


data3RB <- as.data.frame(cbind(fig3RB$Time,fig3RB$OD)); names(data3RB) <- c("time","N")
fit3RB <- fit(data3RB,free=free,fun=log,lower=lower,pch=".",legend=FALSE,tstep=0.1,
              main="red fig 3b lag phase params not ignored",
              ymin = 0, ymax = 1)
summary(fit3RB)

data3GB <- as.data.frame(cbind(fig3GB$Time,fig3GB$OD)); names(data3GB) <- c("time","N")
fit3GB <- fit(data3GB,free=free,fun=log,lower=lower,pch=".",legend=FALSE,tstep=0.1,
              main="green fig 3b lag phase params not ignored",
              ymin = 0, ymax = 1)
summary(fit3GB)

#Get for both bacteria that you cannot get covariance and that the system is singular.
#See : https://stackoverflow.com/questions/50928796/system-is-computationally-singular-reciprocal-condition-number-in-r 
#https://math.stackexchange.com/questions/889425/what-does-determinant-of-covariance-matrix-give 
#https://stackoverflow.com/questions/34700685/calculate-the-inverse-of-a-matrix-system-is-computationally-singular-error 

#Fit is 'good' but parameters show that q0 goes extremely high. --> eliminates lag by making that fraction
#effectively zero.



#######################################################################
#Question 3A: fitting figure 3 A while only having 3 parameters differ#
#######################################################################

#Adapt the example, remove N, because should be the same.
#Use the monoculture model (so model = model), which is defined at the top of the file.

#Let them play with different combinations and think about what that means.
#Do (K, m, v), (K, q0), K alone, (K, v) etc.

s <- c(N=0.124)
p <- c(K=0.6,r=0.4,m=2,q0=0.005,v=2)
free <- c("N",names(p))
differ <- c("K", "v", "m")
totfree <- c(free[!(free %in% differ)], differ, differ); npar <- length(totfree)
cat("Number of free parameters",npar)
lower <- rep(0,npar); lower[which(totfree == "v")] <- 1; lower  # set lower bounds
fitq3aKVMDiffer <- fit(data=list(data3RA,data3GA),free=free,differ=differ,fun=log, odes = model,
             lower=lower,pch=".",legend=FALSE,tstep=0.1,main="red & green K, v and m differing",add=TRUE,
             ymin = 0, ymax = 1)
summary(fitq3aKVMDiffer)
fitq3aKVMDiffer$ssr

#So, when fitting 3A while letting only K, v and m differ
#3R (exptA) params:  K = 0.65; v = 1;   m = 2.48
#3G (exptA) params:  K = 0.53; v = 1.4; m = 0.77 
#qo and r and N (shared) : N = 0.12; q0 = 0.03; r = 0.61

#compare with original:
summary(fit3RA)
fit3RA$ssr
#3R (exptA) : N = 0.13; K = 0.65; v = 1; m = 3.21; q0 = 0.01; r = 0.61

summary(fit3GA)
fit3GA$ssr
#3G (exptA) : N = 0.12; K = 0.52; v = 3.1; m = 0.95; q0 = 0.03; r = 0.34

fit3RA$ssr + fit3GA$ssr
#fit of the model with 2 datasets is only very minimally worse in ssr than the two separate
abs(fitq3aKVMDiffer$ssr-(fit3RA$ssr + fit3GA$ssr))


#Do this also with only K and q0 differing:

differ <- c("K", "q0")
totfree <- c(free[!(free %in% differ)], differ, differ); npar <- length(totfree)
cat("Number of free parameters",npar)
lower <- rep(0,npar); lower[which(totfree == "v")] <- 1; lower  # set lower bounds
fitq3aKQ0Differ <- fit(data=list(data3RA,data3GA),free=free,differ=differ,fun=log, odes = model,
             lower=lower,pch=".",legend=FALSE,tstep=0.1,main="red & green K and q0 differing",add=TRUE,
             ymin = 0, ymax = 1)
summary(fitq3aKQ0Differ)
fitq3aKQ0Differ$ssr

abs(fitq3aKQ0Differ$ssr-(fit3RA$ssr + fit3GA$ssr))
#1.415612 --> so quite some difference if you only allow K and q0 to vary, but still works well enough.

######################################################################################
##Question 3B : fitting red fig 3 A and B together; and green fig 3 A and B together##
##                                (Old question 4a, hence names)                    ##
######################################################################################

fig3RA <- read.csv(paste("Fig3/",exptA,"_R.csv",sep="")) # Red
data3RA <- as.data.frame(cbind(fig3RA$Time,fig3RA$OD)); names(data3RA) <- c("time","N")
fig3GA <- read.csv(paste("Fig3/",exptA,"_G.csv",sep="")) # Green
data3GA <- as.data.frame(cbind(fig3GA$Time,fig3GA$OD)); names(data3GA) <- c("time","N")
fig3RB <- read.csv(paste("Fig3/",exptB,"_R.csv",sep="")) # Red
data3RB <- as.data.frame(cbind(fig3RB$Time,fig3RB$OD)); names(data3RB) <- c("time","N")
fig3GB <- read.csv(paste("Fig3/",exptB,"_G.csv",sep="")) # Green
data3GB <- as.data.frame(cbind(fig3GB$Time,fig3GB$OD)); names(data3GB) <- c("time","N")
fig3RC  <- read.csv(paste("Fig3/",exptC,"_R.csv",sep="")) # Red
data3RC <- as.data.frame(cbind(fig3RC$Time,fig3RC$OD)); names(data3RC) <- c("time","N")
fig3GC  <- read.csv(paste("Fig3/",exptC,"_G.csv",sep="")) # Red
data3GC <- as.data.frame(cbind(fig3GC$Time,fig3GC$OD)); names(data3GC) <- c("time","N")

#A1 and B1 = data3RA data3RB

model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    a <- q0/(q0+exp(-m*t))
    dN <- r*a*N*(1-(N/K)^v)
    return(list(dN))  
  }) 
}


sQ4a <- c(N=0.124) #note that now this estimate is only good for 1 of the 2, the other starts more around ~0.22
pQ4a <- c(K=0.6,r=0.4,m=2,q0=0.005,v=2)
freeQ4a <- c("N",names(pQ4a))
differQ4a <- c("N", "q0", "m")
totfreeQ4a <- c(freeQ4a[!(freeQ4a %in% differQ4a)], differQ4a, differQ4a); nparQ4a <- length(totfreeQ4a)
cat("Number of free parameters",nparQ4a)
lowerQ4a <- rep(0,nparQ4a); lowerQ4a[which(totfreeQ4a == "v")] <- 1; lowerQ4a  # set lower bounds
fitQ4aRed <- fit(data=list(data3RA, data3RB), state = sQ4a, parms = pQ4a, free=freeQ4a,differ=differQ4a,fun=log,
              lower=lowerQ4a, odes = model, pch=".",legend=FALSE,tstep=0.1,main="red fig3A & red fig 3B",add=TRUE,
              ymin = 0, ymax = 1)
summary(fitQ4aRed)
fitQ4aRed$ssr
fitQ4aRed$par

####
####RED
####
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
#   K  6.472e-01  1.618e-03 400.064  < 2e-16 ***
#   r  6.099e-01  3.474e-02  17.558  < 2e-16 ***
#   v  1.000e+00  8.102e-02  12.343  < 2e-16 ***
#   N  1.256e-01  4.554e-04 275.692  < 2e-16 ***
#   q0 1.170e-02  3.303e-03   3.543 0.000408 ***
#   m  3.169e+00  2.381e-01  13.312  < 2e-16 ***
#   N  2.301e-01  1.270e-03 181.261  < 2e-16 ***
#   q0 8.056e+01  1.360e+02   0.592 0.553765    
#   m  5.359e-05  8.244e-01   0.000 0.999948 

#RESIDUAL SUM OF SQUARES: 1.201063


fitQ4aGreen <- fit(data=list(data3GA, data3GB), state = sQ4a, parms = pQ4a, free=freeQ4a,differ=differQ4a,fun=log,
                 lower=lowerQ4a, odes = model, pch=".",
                 legend=FALSE,tstep=0.1,main="green fig3A & green fig 3B",add=TRUE,
                 ymin = 0, ymax = 1)
summary(fitQ4aGreen)
fitQ4aGreen$ssr
fitQ4aGreen$par

####
####Green
####


# Parameters:
#   Estimate Std. Error t value Pr(>|t|)
# K  6.034e-01         NA      NA       NA
# r  5.484e-01         NA      NA       NA
# v  1.000e+00         NA      NA       NA
# N  1.237e-01         NA      NA       NA
# q0 3.195e-02         NA      NA       NA
# m  8.994e-01         NA      NA       NA
# N  2.967e-01         NA      NA       NA
# q0 1.033e+18         NA      NA       NA
# m  4.840e-01         NA      NA       NA


#Seems to work quite well, given the residual sum of squares.
#Best parameters would be the one reported here: any other differences are probably not rooted in reality, they are the same strain from the same stock!







################################################
#Question 4a : competition model per experiment#
################################################


# The 2D model is defined:

model2 <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    a1 <- q01/(q01+exp(-m1*t)); a2 <- q02/(q02+exp(-m2*t))
    dN1 <- r1*a1*N1*(1-(N1/K1)^v1-c2*(N2^v2)/(K1^v1))
    dN2 <- r2*a2*N2*(1-(N2/K2)^v2-c1*(N1^v1)/(K2^v2))
    return(list(c(dN1,dN2)))  
  }) 
}




#get c1 and c2 values per experiment
#Loops over data list about figure 3 to get monoculture parameters
#Then fits model2 with those parameters given, estimating only c1 and c2.
listFitsExptsFig4 = list()

for (i in seq_along(expSeq)) {
  
  print(i); print(expSeq[i])
  
  
  # Retrieve parameters from the monoculture fits of figure 3 list and rename them for the 2D model
  fig3RFitThisExperiment <- listFitsFigure3[[paste0("Fit", names(expSeq)[i], "Fig3R" )]]
  fig3GFitThisExperiment <- listFitsFigure3[[paste0("Fit", names(expSeq)[i], "Fig3G" )]]
  
  pR <- fig3RFitThisExperiment$par[2:length(fig3RFitThisExperiment$par)]
  names(pR) <- paste(names(pR),"1",sep="")
  
  pG <- fig3GFitThisExperiment$par[2:length(fig3GFitThisExperiment$par)]
  names(pG) <- paste(names(pG),"2",sep="")
  
  
  fig4 <- read.csv(paste("Fig4/",expSeq[i],"_RG.csv",sep=""))
  
  data4 <- as.data.frame(cbind(fig4$Time,fig4$OD)); names(data4) <- c("time","OD")
  p4 <- c(pR,pG,c1=1,c2=1)
  initialOD <- data4[1,2]
  s4 <- c(N1=initialOD/2,N2=initialOD/2);s4  # assume the expt was started equally
  free4 <- c("c1","c2")
  fit4 <- fit(data4,odes=model2,tweak="nsol$OD=nsol$N1+nsol$N2",free=free4, parms = p4, state = s4,
              fun=log,lower=0,upper=2,pch=".",legend=TRUE,tstep=0.1,
              main=paste0("Fig 4 fit using monoculture param estimates from fig3R and G for ", names(expSeq[i])),
              ymin = 0, ymax = 1)
  summary(fit4)
  
  listFitsExptsFig4[[paste0("Data", names(expSeq)[i], "Fig4")]] <- data4
  listFitsExptsFig4[[paste0("Fit", names(expSeq)[i], "Fig4")]]  <- fit4
  
}

summary(listFitsExptsFig4$FitexptAFig4)
summary(listFitsExptsFig4$FitexptBFig4)
summary(listFitsExptsFig4$FitexptCFig4)

#Answer:
# expt A -> C1 = 0.59970;   C2 = 1.99984
# expt B -> C1 = 0.0004199749; C2 = 1.6586690075
# expt C -> C1 = 0.7874716192; C2 = 0.0005538458

#What this means:
#ExptA : the Red strain (c1) experiences ~4 times less competition from green (c2) than vice versa (i.e. red is more hindered by its own population than by green; while green is hindered about 2 as much by red as by itself)
#However, in reality: first strain (red) grows much harder than second strain (green).
#c2 not identifiable because green is basically gone by the time competition starts being
#a real pain (near carrying capacity), so the fit just guesses something arbitrarily





###########################################################################################
#Question 4a cont. How much fidelity do you lose in the fit when you simply set c1=c2 = 1?#
###########################################################################################


#to get at SSR or RSS (residual sum of squares), simply take fit$ssr

#change model

model2EqualC <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    a1 <- q01/(q01+exp(-m1*t)); a2 <- q02/(q02+exp(-m2*t))
    dN1 <- r1*a1*N1*(1-(N1/K1)^v1-c*(N2^v2)/(K1^v1))
    dN2 <- r2*a2*N2*(1-(N2/K2)^v2-c*(N1^v1)/(K2^v2))
    return(list(c(dN1,dN2)))  
  }) 
}

#change pars

listFitsExptsFig4EqualC = list()

#rerun with c instead of c1 and c2
for (i in seq_along(expSeq)) {
  
  print(i); print(expSeq[i])
  
  
  fig3RFitThisExperiment <- listFitsFigure3[[paste0("Fit", names(expSeq)[i], "Fig3R" )]]
  fig3GFitThisExperiment <- listFitsFigure3[[paste0("Fit", names(expSeq)[i], "Fig3G" )]]
  
  #Red is N1
  pREqualC <- fig3RFitThisExperiment$par[2:length(fig3RFitThisExperiment$par)]
  names(pREqualC) <- paste(names(pREqualC),"1",sep="")
  
  #Green is N2
  pGEqualC <- fig3GFitThisExperiment$par[2:length(fig3GFitThisExperiment$par)]
  names(pGEqualC) <- paste(names(pGEqualC),"2",sep="")
  
  
  fig4 <- read.csv(paste("Fig4/",expSeq[i],"_RG.csv",sep=""))
  
  data4 <- as.data.frame(cbind(fig4$Time,fig4$OD)); names(data4) <- c("time","OD")
  p <- c(pREqualC,pGEqualC,c=1)
  initialOD <- data4[1,2]
  s <- c(N1=initialOD/2,N2=initialOD/2);s  # assume the expt was started equally
  free <- c("c")
  fit4 <- fit(data4,odes=model2EqualC,tweak="nsol$OD=nsol$N1+nsol$N2",free=free,
              fun=log,lower=0,upper=2,pch=".",legend=TRUE,tstep=0.1,
              main=paste0("Fig 4 fitting with monoculture parameters and equal c for ",names(expSeq)[i] ),
              ymin = 0, ymax = 1)
  summary(fit4)
  
  listFitsExptsFig4EqualC[[paste0("Data", names(expSeq)[i], "Fig4")]] <- data4
  listFitsExptsFig4EqualC[[paste0("Fit", names(expSeq)[i], "Fig4")]]  <- fit4
  
}

#see new c parameter
summary(listFitsExptsFig4EqualC$FitexptAFig4)
summary(listFitsExptsFig4EqualC$FitexptBFig4)
summary(listFitsExptsFig4EqualC$FitexptCFig4)

#ExptA fig 4 c ==> 0.60
#ExptB fig 4 c ==> 0.97
#ExptC fig 4 c ==> 0.79


#calculate differences in residual sum of squares:

ssrWithDiffC = c(listFitsExptsFig4$FitexptAFig4$ssr, listFitsExptsFig4$FitexptBFig4$ssr,
                 listFitsExptsFig4$FitexptCFig4$ssr)
ssrEqualC    = c(listFitsExptsFig4EqualC$FitexptAFig4$ssr, listFitsExptsFig4EqualC$FitexptBFig4$ssr,
                 listFitsExptsFig4EqualC$FitexptCFig4$ssr)
dataFrameSSR = data.frame(cbind(ssrWithDiffC, ssrEqualC), row.names = c("exptA", "exptB", "exptC"))
dataFrameSSR[,"AbsDifferenceSSR"] = abs(dataFrameSSR[,"ssrWithDiffC"]- dataFrameSSR[,"ssrEqualC"])
print(dataFrameSSR)

#Negligible differences in all experiments.




######################################################################
#Unused question 1. Old question 3c. Fit figure 3 and 4 data together#
######################################################################

#Want to fit the full 2D model
#You have shown that differing only m and K is fine, so you want them to differ in those 2.

model2 <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    a1 <- q01/(q01+exp(-m1*t)); a2 <- q02/(q02+exp(-m2*t))
    v1 <- max(v1, 1); v2 <- max(v2, 1)
    dN1 <- r1*a1*N1*(1-(N1/K1)^v1-c*(N2^v2)/(K1^v1))
    dN2 <- r2*a2*N2*(1-(N2/K2)^v2-c*(N1^v1)/(K2^v2))

    return(list(c(dN1,dN2)))  
  }) 
}

# Retrieve parameters from the "fit$par" list and rename them for the 2D model

pR <- fit3RA$par[2:length(fit3RA$par)]; names(pR) <- paste(names(pR),"1",sep="")
pG <- fit3GA$par[2:length(fit3GA$par)]; names(pG) <- paste(names(pG),"2",sep="")


fig4DataExptA = listFitsExptsFig4$DataexptAFig4
initialNFig4ExptAPerStrain = fig4DataExptA$OD[1]/2

#need to rename N to N1 and N2 so the fit understands what is what. Red = N1; Green = N2.
head(data3RA)
colnames(data3RA) <- c("time", "N1")
head(data3GA)
colnames(data3GA) <- c("time", "N2")
head(fig4DataExptA)

pQuestion3C <- c(pR, pG, c=1)
freeQuestion3C <- c(names(pQuestion3C))

differQuestion3C = c("K1", "K2", "m1", "m2")
#should fix the N1 and N2. Red is not present in the green monoculture and vice versa
fixedQuestion3C = list(N1 = c(0.124, 0, initialNFig4ExptAPerStrain),
                       N2 = c(0, 0.124, initialNFig4ExptAPerStrain)
                       )
#Not used, since it takes this from fixed, but need to set it so it know that these are the state variables
stateQuestion3C <- c(N1 = 0, N2 = 0)
#fit 3 experiments together, so you get 1 time the parameters fitted over all, and 3 times the parameters
#fitted separately per dataset
totfreeQuestion3C <- c(freeQuestion3C[!(freeQuestion3C %in% differQuestion3C)],
                       differQuestion3C,
                       differQuestion3C,
                       differQuestion3C);

#the lower bound should be 1 for both v parameters.
nparQuestion3C = length(totfreeQuestion3C)
lowerQuestion3C <- rep(0,nparQuestion3C);
lowerQuestion3C[which(totfreeQuestion3C == "v1")] <- 1
lowerQuestion3C[which(totfreeQuestion3C == "v2")] <- 1

#Usage of Nelder-Mead makes sure no errors appear in integration, don't know why
fitQuestion3C <- fit(data= list(data3RA,data3GA, fig4DataExptA), parms = pQuestion3C,
                     free = freeQuestion3C, differ = differQuestion3C, state = stateQuestion3C,
                     fixed = fixedQuestion3C, fun = log, odes = model2, lower = lowerQuestion3C,
                     method = "Nelder-Mead",
                     tweak="nsol$OD=nsol$N1+nsol$N2",
                     pch=".",legend=TRUE,tstep=0.1,
                     main="Fitting red and green figure 3A and figure 4a all together (letting only K1 and K2 and m1 and m2 vary)",
                     add=TRUE, ymin = 0, ymax = 1)

summary(fitQuestion3C)




##############################################################################################
##Unused Question 2, old question 4b. Do you expect c1 and c2 to differ between fig4A and B?##
##############################################################################################

#Simple answer is no: should be the same strains, in the same experimental conditions. 
#Test by fitting 4A and 4B together, allowing only c1 and c2 (or even c) to differ!

#Best possible way to do this is to:
#fit the data from 3ARed and 3BRed together and save those parameters
#and fit the data from 3AGreen and 3BGreen together and save those parameters
#-->BOTH OF THESE ARE DONE IN QUESTION 3B!
#then fit fig4 A and B together using the monoculture parameters estimated from those two, using
#the lag phase q0 and m values (because you need those?)


#bedoeling: fit 3a, krijg 14 parameterwaardes, fit daarmee figuur 4a data. Haal uit monocultuur alle
#groeiparameters, 


#Now fit the data of fig4A and B together. Make sure that the right q0 and m are used per 
#experiment (in A, no lag phase or a very large q0 that eliminates it; in B the other)


model2 <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    a1 <- q01/(q01+exp(-m1*t)); a2 <- q02/(q02+exp(-m2*t))
    dN1 <- r1*a1*N1*(1-(N1/K1)^v1-c2*(N2^v2)/(K1^v1))
    dN2 <- r2*a2*N2*(1-(N2/K2)^v2-c1*(N1^v1)/(K2^v2))
    return(list(c(dN1,dN2)))  
  }) 
}

#give the parameters names. Taken from answer of question 3c above.
#Again: 3c was once 4a, hence the names on the right here.
parametersGreenFig3AandB <- fitQ4aGreen$par
parametersRedFig3AandB   <- fitQ4aRed$par

#first 3 parameters were fitted shared, last 6 were fitted separately
names(parametersRedFig3AandB) <- paste0(names(parametersRedFig3AandB), "1")

#rename those fitted in exptA with _A and those in B with _B, for easy setting in fixed list below.
names(parametersRedFig3AandB)[4:6] <- paste0(names(parametersRedFig3AandB)[4:6], "_A")
names(parametersRedFig3AandB)[7:9] <- paste0(names(parametersRedFig3AandB)[7:9], "_B")

names(parametersGreenFig3AandB) <- paste0(names(parametersGreenFig3AandB), "2")
names(parametersGreenFig3AandB)[4:6] <- paste0(names(parametersGreenFig3AandB)[4:6], "_A")
names(parametersGreenFig3AandB)[7:9] <- paste0(names(parametersGreenFig3AandB)[7:9], "_B")

parametersGreenFig3AandB; parametersRedFig3AandB
#green = 2, red = 1

#pR and pG are the names for the parameters of green and red estimated from their respective
#monoculture fits in experiment A
pQuestion4b <- c(pR, pG, c1=1, c2 = 1)
freeQuestion4b <- c("c1", "c2")


#right now, I took the q and m values for each strain that correspond to the fit when there is a lag phase, i.e. in expt A.
fixedQuestion4b = list(q02 = c(parametersGreenFig3AandB["q02_A"], parametersGreenFig3AandB["q02_B"]), 
                       m2  = c(parametersGreenFig3AandB["m2_A"], parametersGreenFig3AandB["m2_B"]),
                       N2  = c(0.06, 0.1281),
                       q01 = c(parametersRedFig3AandB["q01_A"], parametersRedFig3AandB["q01_B"]), 
                       m1  = c(parametersRedFig3AandB["m1_A"], parametersRedFig3AandB["m1_B"]),
                       N1  = c(0.06, 0.1281)
                       )
stateQuestion4b <- c(N1 = 0, N2 = 0)

fitQuestion4b <- fit(data= list(listFitsExptsFig4$DataexptAFig4,
                                listFitsExptsFig4$DataexptBFig4),
                     parms = pQuestion4b,
                     free = freeQuestion4b, state = stateQuestion4b,
                     fixed = fixedQuestion4b, fun = log, odes = model2,
                     tweak="nsol$OD=nsol$N1+nsol$N2", lower = 0, method = "Nelder-Mead",
                     pch=".",legend=TRUE,tstep=0.1,main="Fitting red and green of figure 4 a and b together (using monoculture params from 3A and B)",
                     add=TRUE, ymin = 0, ymax = 1)

summary(fitQuestion4b)
#I had for expt A and B:
# expt A -> C1 = 0.59970;   C2 = 1.99984
# expt B -> C1 = 3.278e-01; C2 = 5.368e-04






###############################################
#Unused question 3, old question 5. NOT DONE  # 
###############################################

#Do something with figure 5?


p["c1"] <- fit4$par["c1"]; p["c2"] <- fit4$par["c2"]; p # Retrieve parameters from "fit$par"
initialF <- fig5R$freq_mean[1]
s <- c(N1=initialF*initialOD,N2=(1-initialF)*initialOD)
nsol <- run(6,0.1,odes=model2,tweak="nsol$fR=nsol$N1/(nsol$N1+nsol$N2);nsol$fG=1-nsol$fR",table=TRUE)
plot(nsol$time,nsol$fR,ylim=c(0,1),type="l",col="red")
points(fig5R$time, fig5R$freq_mean, col="red")
lines(nsol$time,nsol$fG,col="darkGreen")
points(fig5G$time, fig5G$freq_mean, col="darkGreen")

# Here the paper ends



##############################################
#Unused question 4, old question 6. NOT DONE #
##############################################


# This is an example of fitting the population genetics model to frequency data

model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    df <- r*s*f*(1-f)
    return(list(df))  
  }) 
}

s <- c(f=0.5); p <- c(r=1,s=0.1); free=c("f","s")

data <- as.data.frame(cbind(fig5R$time,fig5R$freq_mean)); names(data) <- c("time","f")
fitq4d1 <- fit(data,free=free)
data <- as.data.frame(cbind(nsol$time,nsol$fR)); names(data) <- c("time","f")
fitq4d2 <- fit(data,free=free)




###################################################
#Example code for fitting data sets simultaneously#
###################################################


s <- c(N=0.124)
p <- c(K=0.6,r=0.4,m=2,q0=0.005,v=2)
free <- c("N",names(p))
differ <- c("N","K","v","m")
totfree <- c(free[!(free %in% differ)], differ, differ); npar <- length(totfree)
cat("Number of free parameters",npar)
lower <- rep(0,npar); lower[which(totfree == "v")] <- 1; lower  # set lower bounds
fitq1 <- fit(data=list(data3RA,data3GA),free=free,differ=differ,fun=log,lower=lower,pch=".",legend=FALSE,tstep=0.1,main="red & green",add=TRUE)
summary(fitq1)

