# The 3 experiments are indexed by their date:
expt <- "2015-11-18"  #Experiment A
#expt <- "2015-12-14"  #Experiment B
#expt <- "2016-01-06"  #Experiment C

# First read and plot all data:

fig3R <- read.csv(paste("Fig3/",expt,"_R.csv",sep="")) # Red
fig3G <- read.csv(paste("Fig3/",expt,"_G.csv",sep="")) # Green
plot(fig3R$Time, fig3R$OD, ylim=c(0,0.8),col="red",pch=".",xlab="Time (hr)",ylab="OD")
points(fig3G$Time, fig3G$OD, ylim=c(0,0.8),col="green",pch=".")

fig4 <- read.csv(paste("Fig4/",expt,"_RG.csv",sep=""))
plot(fig4$Time, fig4$OD, ylim=c(0,0.8), col="blue", pch=".", xlab="Time (hr)", ylab="Total OD")

fig5 <- read.csv(paste("Fig5/flow_df_",expt,".csv",sep=""))
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

s <- c(N=0.124)
p <- c(K=0.6,r=0.4,m=2,q0=0.005,v=2)
free <- c("N",names(p))
lower <- rep(0, length(free)); lower[match("v",free)] <- 1; lower  # set lower bounds
data3R <- as.data.frame(cbind(fig3R$Time,fig3R$OD)); names(data3R) <- c("time","N")
fit3R <- fit(data3R,free=free,fun=log,lower=lower,pch=".",legend=FALSE,tstep=0.1,main="red")
summary(fit3R)

data3G <- as.data.frame(cbind(fig3G$Time,fig3G$OD)); names(data3G) <- c("time","N")
fit3G <- fit(data3G,free=free,fun=log,lower=lower,pch=".",legend=FALSE,tstep=0.1,main="green")
summary(fit3G)

# The 2D model is defined:

model2 <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    a1 <- q01/(q01+exp(-m1*t)); a2 <- q02/(q02+exp(-m2*t))
    dN1 <- r1*a1*N1*(1-(N1/K1)^v1-c2*(N2^v2)/(K1^v1))
    dN2 <- r2*a2*N2*(1-(N2/K2)^v2-c1*(N1^v1)/(K2^v2))
    return(list(c(dN1,dN2)))  
  }) 
}

# Retrieve parameters from the "fit$par" list and rename them for the 2D model

pR <- fit3R$par[2:length(fit3R$par)]; names(pR) <- paste(names(pR),"1",sep="")
pG <- fit3G$par[2:length(fit3G$par)]; names(pG) <- paste(names(pG),"2",sep="")

data4 <- as.data.frame(cbind(fig4$Time,fig4$OD)); names(data4) <- c("time","OD")
p <- c(pR,pG,c1=1,c2=1); p["v1"] <- max(1,p["v1"]); p["v2"] <- max(1,p["v2"]); p
initialOD <- data4[1,2]
s <- c(N1=initialOD/2,N2=initialOD/2);s  # assume the expt was started equally
free <- c("c1","c2")
fit4 <- fit(data4,odes=model2,tweak="nsol$OD=nsol$N1+nsol$N2",free=free,fun=log,lower=0,upper=2,pch=".",legend=FALSE,tstep=0.1,main="blue")
summary(fit4)

p["c1"] <- fit4$par["c1"]; p["c2"] <- fit4$par["c2"]; p # Retrieve parameters from "fit$par"
initialF <- fig5R$freq_mean[1]
s <- c(N1=initialF*initialOD,N2=(1-initialF)*initialOD)
nsol <- run(6,0.1,odes=model2,tweak="nsol$fR=nsol$N1/(nsol$N1+nsol$N2);nsol$fG=1-nsol$fR",table=TRUE)
plot(nsol$time,nsol$fR,ylim=c(0,1),type="l",col="red")
points(fig5R$time, fig5R$freq_mean, col="red")
lines(nsol$time,nsol$fG,col="darkGreen")
points(fig5G$time, fig5G$freq_mean, col="darkGreen")

# Here the paper ends


# This is an example showing how to fit two data sets simultaneously

s <- c(N=0.124)
p <- c(K=0.6,r=0.4,m=2,q0=0.005,v=2)
free <- c("N",names(p))
differ <- c("K","v","m")
totfree <- c(free[!(free %in% differ)], differ, differ); npar <- length(totfree)
cat("Number of free parameters",npar)
lower <- rep(0,npar); lower[which(totfree == "v")] <- 1; lower  # set lower bounds
fitq1 <- fit(data=list(data3R,data3G),free=free,differ=differ,fun=log,lower=lower,pch=".",legend=FALSE,tstep=0.1,main="red & green",add=TRUE)
summary(fitq1)

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
