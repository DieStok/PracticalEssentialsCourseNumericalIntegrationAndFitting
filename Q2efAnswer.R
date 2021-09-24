#### Essential skills tutorial: Bacteria growing on a resource ####

### 1. ODE solving using different methods ###
if(!require(pacman)) {
  install.packages("pacman"); require(pacman)}
p_load(tidyverse, cowplot)
library(deSolve)


# Define the model in format acceptable for deSolve (and Grind)
model <- function(t, state, parms) {
  # state <- ifelse(state<0,0,state)
  with(as.list(c(state,parms)), {
    f <- v*R/(k+R)
    dtR <- - f*e*B
    dtB <- f*B
    return(list(c(dtR, dtB)))  
  }) 
} 

#If you uncomment # state <- ifelse(state<0,0,state), then
#if you get negative values because of overshooting 0,
#you correct that mistake before doing the next 
#numerical integration step

# Parameters
p <- c(e=5e-7,v=1.4,k=1)
# Initial conditions
#Note: This HAS to be in the order in which you define the equations in the list of functions
s <- c(R=350,B=1e3)

# Time vector
timestep <- 0.5     
time <- seq(0,12,by=timestep)

# Solve the system using the Euler method
startTiming = proc.time()
out <- ode(y=s, times=time, func=model, parms=p, method="euler", hini=timestep)
out_euler <- data.frame(out)
timeEuler = proc.time() - startTiming

# Solve the system using Runge-Kutta 23, with dynamical time step
startTiming = proc.time()
out <- ode(y=s, times=time, func=model, parms=p, method="ode23", hini=timestep)
out_rk23 <- data.frame(out)
timeRKMethod = proc.time() - startTiming

# Plot results base R
dev.off()
par(mfrow=c(2,1))
plot(time,out_euler$R, type='l', col='red', lwd=3)
lines(time,out_rk23$R, type='l', col='orange', lwd=3)
plot(time,out_euler$B, type='l', col='blue', lwd=3, log="y")  # Note: logarithmic y-axis
lines(time, out_rk23$B, type='l', col='cyan', lwd=3)

#Plot results ggplot --> arranged with cowplot library (see here: https://wilkelab.org/cowplot/index.html)
dataFramePlotsEulerVsRK = dplyr::bind_rows(out_euler, out_rk23) %>%
  mutate(method = c(rep("Euler", nrow(out_euler)),
                    rep("Runge-Kutta23", nrow(out_rk23)))
         )

Rplot = dataFramePlotsEulerVsRK %>%
ggplot(aes(x = time, y = R, colour = method)) +
  geom_point() + 
  theme_bw() +
  # we set the left and right margins to 0 to remove 
  # unnecessary spacing in the final plot arrangement.
  theme(plot.margin = margin(6, 0, 6, 0))

Bplot = dataFramePlotsEulerVsRK %>%
  ggplot(aes(x = time, y = B, colour = method)) +
  geom_point() + 
  theme_bw() +
  scale_y_continuous(trans = "log",
                     breaks = c(2, 20, 200, 2000, 20000, 200000, 2000000, 20000000, 200000000, 2000000000)
                     ) +
  theme(plot.margin = margin(6, 0, 6, 0))

legend = get_legend(
  # create some space to the left of the legend
  Rplot + theme(legend.box.margin = margin(0, 0, 0, 12))
)

rowPlots = cowplot::plot_grid(Rplot + theme(legend.position = "none"),
                       Bplot + theme(legend.position = "none"),
                       labels = c("A", "B"),
                       align = "vh",
                       hjust = -1,
                       nrow = 2)

plot_grid(rowPlots, legend, rel_widths = c(3, 1))

#what about calculation time?
print(paste0("Calculation Euler for timestep: ", timestep))
print(timeEuler)
print(paste0("Calculation Runge-Kutta23 for timestep: ", timestep))
print(timeRKMethod)


