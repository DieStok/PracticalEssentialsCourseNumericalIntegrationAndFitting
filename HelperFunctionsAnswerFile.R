#dependencies
if(!require(pacman)) {
  install.packages("pacman"); require(pacman)}
p_load(tidyverse, microbenchmark, RColorBrewer)


#functionalising the calculation loop.

runSimpleNumericalIntegration = function(timeStep = 0.5, maxTime = 100, performMicrobenchmark = FALSE) {


#for time step is 5 the plot looks as follows: 
timestep      = timeStep
time          = seq(0, maxTime, by = timestep)

# initial condition
x0 <- 0.1
## The function to be integrated (right-hand expression of the derivative above)
f <- function(x){x * (1.-x)}

## An empty R vector to store the results
x <- c()
## Store the initial condition in the first position of the vector
x[1] <- x0

# loop over time: approximate the function at each time step
ptm <- proc.time()

for (i in 1:(length(time)-1)){
  x[i+1] = x[i] + timestep * f(x[i])
}

timeIntegration = proc.time() - ptm
if (performMicrobenchmark == TRUE) {
  
  performEuler = function(banana = x, time, timestep) {
    
    for (i in 1:(length(time)-1)){
      banana[i+1] = banana[i] + timestep * f(banana[i])
    }
    banana
  }
  
  
  
  outputBenchmark = microbenchmark(performEuler, times = 10, unit = "s")
  
} else {
  
  outputBenchmark = NA
  
}

xAnalytic = 0.1 * exp(time)/(1+0.1*(exp(time)-1.))


#print(time)
#print(x)
#print(xAnalytic)
dataFramePlot <- data.frame(x_value = c(x, xAnalytic),
                                time = time,
                                solution_method = c(rep("Approximation", length(x)), rep("Analytical", length(x)))
)

plot = ggplot2::ggplot(data = dataFramePlot,
                aes (x = time, y = x_value, colour = solution_method)) + 
  geom_point(data = dataFramePlot %>% subset(solution_method == "Approximation"), size = 2) +
  geom_line(data = dataFramePlot %>% subset(solution_method != "Approximation"), size = 1.2) +
  theme_bw() +
  ggtitle(paste0("Approximate versus analytical solution of f(x) = x*(1-x) for time step of ", timestep, "")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_brewer(palette = "Paired")

show(plot)

#This is reflected in the sum of differences:
differences      = abs(x - 0.1 * exp(time)/(1+0.1*(exp(time)-1.)))
totalDifferences = sum(differences)
print(paste0("Total differences between numerically integrated and analytic solution if timestep is ",
             timestep, ": ",
             totalDifferences)
      )

return(list(data = dataFramePlot,
            plot = plot,
            timeExpended    = timeIntegration,
            Totdifferences  = totalDifferences,
            benchMarkData   = outputBenchmark
            )
       )

}

runEulerBacteriaResource      = function(timeStep = 0.5, maxTime = 100, performMicroBenchMark = FALSE, plotType = "point") {
  
  time     = seq(0, maxTime, by = timeStep)
  
  R0 = 350
  B0 = 1000
  f  = function(R,B){-e*(v*R/(k+R))*B}
  g  = function(R,B){(v*R/(k+R)*B)}
  v  = 1.4
  e  = 0.0000007
  k  = 1
  
  ## An empty vector to store the results
  R <-c()
  B <-c()
  ## Store the initial condition in the first position of the vector
  R[1] <- R0
  B[1] <- B0
  # loop over time: approximate the function at each time step
  startTime = proc.time()
  for (i in 1:(length(time)-1)){
    R[i+1]=R[i]+timeStep*f(R[i],B[i])
    B[i+1]=B[i]+timeStep*g(R[i],B[i])
  }
  calculationTime = proc.time() - startTime
  
  #ggplots
  dataFramePlotsRB = data.frame(x = time, R = R, B = B)
  
  if (plotType == "point") {
  #R
  plotR = dataFramePlotsRB %>%
    ggplot(aes(x = x, y = R)) + geom_point() +
    theme_bw() +
    ggtitle(paste0("Approximate solution of R for time step of ", timeStep, "")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  #B
  plotB = dataFramePlotsRB %>%
    ggplot(aes(x = x, y = B)) + geom_point() +
    theme_bw() +
    ggtitle(paste0("Approximate solution of B for time step of ", timeStep, "")) +
    theme(plot.title = element_text(hjust = 0.5))
  } else if (plotType == "line") {
    
    #R
    plotR = dataFramePlotsRB %>%
      ggplot(aes(x = x, y = R)) + geom_line() +
      theme_bw() +
      ggtitle(paste0("Approximate solution of R for time step of ", timeStep, "")) +
      theme(plot.title = element_text(hjust = 0.5))
    
    #B
    plotB = dataFramePlotsRB %>%
      ggplot(aes(x = x, y = B)) + geom_line() +
      theme_bw() +
      ggtitle(paste0("Approximate solution of B for time step of ", timeStep, "")) +
      theme(plot.title = element_text(hjust = 0.5))
    
    
  } else {
    stop("Error, plotType should be one of line or point")
  }
  
  
  
  show(plotR)
  show(plotB)
  print(paste0("Calculation time for time step ", timeStep,
               " for ", maxTime, " steps total: ", calculationTime[3] 
               )
        )
  
  
  return( list(data = dataFramePlotsRB, plots = list(plotR, plotB)))
  
}

readInDataFigureThreeRamPaper = function() {
  
  
  #model to fit
  model <- function(t, state, parms) {
    with(as.list(c(state,parms)), {
      a <- q0/(q0+exp(-m*t))
      dN <- r*a*N*(1-(N/K)^v)
      return(list(dN))  
    }) 
  }
  
  # The 3 experiments are indexed by their date:
  exptA <- "2015-11-18"  #Experiment A
  exptB <- "2015-12-14"  #Experiment B
  exptC <- "2016-01-06"  #Experiment C
  
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
  return(listFitsFigure3)
}


fitExpt3AWithDiffer           = function(differ = c("K", "v", "m")) {
  
  cat("   \n")
  cat("---\n")
  print("Differ for this call: ")
  print(differ)
  cat("---\n")
  cat("   \n")
  #set up state, parameters, and what to fit separately
  initialState  <- c(N=0.124)
  initialParams <- c(K=0.6,r=0.4,m=2,q0=0.005,v=2)
  free <- c("N",names(initialParams))
  
  #the total number of parameters to fit = 2 * those in differ + all parameters in free that should be fitted together for the two datasets
  totfree <- c(free[!(free %in% differ)], differ, differ)
  npar <- length(totfree)
  cat("Number of free parameters",npar, "\n")
  
  #set the lower bound to 0 for all, then to 1 for all instances of v
  lower <- rep(0,npar); lower[which(totfree == "v")] <- 1; lower  # set lower bounds
  
  
  fitDiffer <- fit(data=list(dataAndFitsFigureThree$DataexptAFig3R, dataAndFitsFigureThree$DataexptAFig3G), 
                      free=free, differ=differ, fun=log, odes = model, state = initialState,
                      parms = initialParams, lower=lower, pch=".",legend=FALSE, tstep=0.1, 
                      main=paste("red & green", paste(differ, collapse = ", "), "differing"), add=TRUE, ymin = 0, ymax = 1)
  summary(fitDiffer)
  fitDiffer$ssr
  

  
  #compare with original two fits separately:
  summary(dataAndFitsFigureThree$FitexptAFig3R)
  dataAndFitsFigureThree$FitexptAFig3R$ssr
  
  summary(dataAndFitsFigureThree$FitexptAFig3G)
  dataAndFitsFigureThree$FitexptAFig3G$ssr
  
  #fit of the model with 2 datasets is only very minimally worse in ssr than the two separate
  #note: should of course add up the residuals of the fit of Red and Green separately for a fair comparison!
  differencesSSR = abs(fitDiffer$ssr-(dataAndFitsFigureThree$FitexptAFig3R$ssr + dataAndFitsFigureThree$FitexptAFig3G$ssr))
  print(paste0("Difference in ssr between fitted together with parameter(s): ", paste(differ, collapse = ", "), " differing : ", differencesSSR))
  
  return(list(fit = fitDiffer, diffSSR = differencesSSR))
}
