# Script used in the research paper "Bayesian approach and extreme value
# theory in economic analysis of forestry projects". Any doubts or
# comments should be addressed to thiagotsalles@gmail.com

# This script requires OpenBUGS to be installed on the computer

# Libraries
library(ggplot2)    # for plots
library(reshape2)   # for melt
library(cowplot)    # for get_legend
library(R2OpenBUGS) # for bugs
library(coda)       # for as.mcmc
library(Kendall)    # for MannKendall

# Function for calculating RMSE
RMSE = function(res) {
  n = length(res)
  sq_err = 0
  for (i in 1:n) {
    sq_err = sq_err + res[i] ^2
  }
  rmse = sqrt(sq_err / n)
}

#### =========================== VOLUME =========================== ####
# DATA
# Prior information for volume
volume.p1 <- read.csv("prior_info1.csv", header=T)
volume.p2 <- read.csv("prior_info2.csv", header=T)

# Sample for volume
volume.sampl <- read.csv("sample.csv", header=T)

# Binding and labeling data
p1_long <- melt(volume.p1, id="age")
p2_long <- melt(volume.p2, id="age")
p1_long$name <- "Prior 1"
p2_long$name <- "Prior 2"
volume.sampl$name <- "Sample"
vol.bind <- rbind(p1_long, p2_long)

# GRAPHS
# Some formatting
formatPlot = function(ggpobj) {
  ggpobj +
  scale_y_continuous(
    breaks=seq(0, 600, 200), minor_breaks=seq(0, 600, 100),
    limits=c(0, 625)) +
  theme_grey(base_size=12) +
  theme(axis.text=element_text(size=12, color="black"),
        plot.title=element_text(face="bold"))
}

# Graph for sample
sample.graph <- ggplot(
  volume.sampl,aes(x=age, y=value, group=variable)) +
  geom_point(size=2) +
  theme(axis.ticks=element_line(color="black")) +
  labs(title="Sample", x="Age (years)", y=bquote('Volume'~(m^3~ha^-1)))

sample.graph <- formatPlot(sample.graph)

# Graph for prior information
prior.graph <- ggplot(
  vol.bind, aes(x=age, y=value, group=variable, shape=name)) +
  geom_point(size=2, fill="white") + 
  geom_point(aes(shape="Sample"), size=2, alpha="0") +  # dummy
  scale_shape_manual(
    breaks=c("Sample", "Prior 1", "Prior 2"),
    values=c(18, 21, 19)) +
  labs(title="Prior information", x="Age (years)", y=element_blank())

prior.graph <- formatPlot(prior.graph) +
  theme(axis.ticks.y=element_blank(), axis.text.y=element_blank())

# Legend for graphs
legend.volume <- get_legend(prior.graph +
                            theme(legend.title=element_blank(),
                                  legend.position="top",
                                  legend.direction="horizontal",
                                  legend.text=element_text(size=12)))

# Joining sample and prior graphs
volume.graph <- plot_grid(sample.graph + theme(legend.position="none"),
                          prior.graph + theme(legend.position="none"),
                          rel_widths=c(1.05, 0.95))

# Show volume graph
plot_grid(volume.graph, legend.volume, nrow=2, rel_heights = c(1, 0.2))

  
# FREQUENTIST FIT FOR PRIOR PARAMETERS
vol.split <- split(vol.bind, vol.bind$variable)
b0.list <- c()
b1.list <- c()
for(i in 1 : length(vol.split)) {
  x = vol.split[[i]]$age
  y = vol.split[[i]]$value
  schumacher = nls(y ~ exp(b0 - b1 / x), start=c(b0=5, b1=5))
  b0.list[i] <- coefficients(schumacher)[[1]]
  b1.list[i] <- coefficients(schumacher)[[2]]
}

mu0.b0 <- mean(b0.list); mu0.b0
tau0.b0 <- 1 / sd(b0.list) ^ 2; tau0.b0
mu0.b1 <- mean(b1.list); mu0.b1
tau0.b1 <- 1 / sd(b1.list) ^ 2; tau0.b1

# BAYESIAN FIT
# Data
x <- volume.sampl$age
y <- volume.sampl$value

# Items required to fit the model with OpenBUGS
ages <- seq(1, 6, 1) # ages for estimating volume by MCMC
n.bayes <- length(x)
data <- list("x", "y", "n.bayes", "ages")
theta <- c("b0", "b1", "tau2")
inits <- function() {list(b0=5, b1=5, tau2=1.0)}

# Creating a text file to run the "bugs" function in OpenBUGS software
sink("schumacher.txt")
cat("
    model{
    
    # PRIOR DISTRIBUTIONS
    b0 ~ dnorm(6.8026, 3.3033)
    b1 ~ dnorm(5.7388, 1.2481)
    tau2 ~ dgamma(0.001, 0.001)
    
    # MODEL
    for(i in 1 : n.bayes) {
    y[i] ~ dnorm(mu[i], tau2)
    mu[i] <- exp(b0 - b1 / x[i])
    }
    
    # ESTIMATES FOR WEIGHT
    for(i in 1 : 6) {
    yp[i] <- exp(b0 - b1 / ages[i])
    }

    # R2B CALCULATION
    sigma2 <- 1 / tau2
    sy2 <- pow(sd(y[]), 2)
    R2B <- 1 - sigma2 / sy2

    }", fill=TRUE)
sink()

# Model fit. On Windows, OpenBUGS will open to run the "bugs" function.
# OpenBUGS has to be closed after it finishes running its algorithm so
# the results can return to R.
schumacher.vol = bugs(data, inits,
                      model.file="schumacher.txt",
                      parameters=c(theta, "yp", "R2B"),
                      n.chains=1, n.iter=12000, n.burnin=5000, n.thin=4,
                      codaPkg=FALSE)

# Convergence tests
schumacher.vol.mcmc <- as.mcmc(schumacher.vol$sims.matrix[, ])
geweke.diag(schumacher.vol.mcmc)
raftery.diag(schumacher.vol.mcmc)
heidel.diag(schumacher.vol.mcmc)

# Results table, with parameters and weight estimates
hpd.interval <- HPDinterval(schumacher.vol.mcmc)
results.schum <- data.frame(
  "Estimate"=schumacher.vol$summary[, 1],
  "HPD lower limit"=hpd.interval[, 1],
  "HPD upper limit"=hpd.interval[, 2])
results.schum <- round(results.schum, 3)

# RMSE
b0 <- schumacher.vol$summary[, 1][1]
b1 <- schumacher.vol$summary[, 1][2]
y.res <- y - exp(b0 - b1 / x)
rmse.schum <- RMSE(y.res)

#### =========================== PRICE ============================ ####
# DATA
prices <- read.csv("price.csv", header=T)
prices <- melt(prices, id="year")
prices.min <- subset(prices, variable == "Minimum")[, 3]
prices.max <- subset(prices, variable == "Maximum")[, 3]

# Tests for independance and trend
Box.test(as.numeric(prices.min), type=("Ljung-Box"))
MannKendall(as.numeric(prices.min))
Box.test(as.numeric(prices.max), type=("Ljung-Box"))
MannKendall(as.numeric(prices.max))

# GRAPH 
price.graph <- ggplot(
  prices, aes(x=year, y=value, linetype=variable)) +
  geom_line() +
  labs(x="Year", y=bquote("Price (R$"~m^3*")")) +
  scale_linetype_manual(values=c(2, 1)) +
  scale_y_continuous(limits=c(37, 62)) +
  scale_x_continuous(breaks=seq(2006, 2017, 2)) +
  theme_grey(base_size=12) +
  theme(axis.text = element_text(size=12, color="black"),
        legend.title=element_blank(),
        legend.position="top",
        legend.text=element_text(size=12))

# Show price graph
price.graph

# PRELIMINAR FIT MIN
# Inputs for minimum
trend.p <- seq(1, 9)
min.p.sample <- prices.min[1:9]
data.p.min <- list(x=trend.p, y=-min.p.sample, n=length(min.p.sample))
inits.p.min <- function(){ list(alpha0=-50, beta=2, alpha1=0.1)}

# Creating a text file to run the "bugs" function in OpenBUGS software
sink("gumbel_p_min.txt")
cat("
    model {
    
    # PRIOR DISTRIBUTIONS
    alpha0 ~ dnorm(0.0, 0.000001)
    beta ~ dnorm(0.0, 0.000001)
    alpha1 ~ dnorm(0.0, 0.0001)
    
    # LIKELIHOOD
    for (i in 1:n) {
    y[i] ~ dgumbel(alpha[i], beta)
    alpha[i] <- alpha0 + alpha1 * x[i]
    }
    
    # PREDICTIVE
    for(i in 2:3) {
    min_p[i] <- ((alpha0 + alpha1 * (9 + i)) - 
                 beta * log(-log(1 - (1 / i))))
    }
    
    }", fill=TRUE)
sink()

# Distribution fit. On Windows, OpenBUGS will open to run the "bugs"
# function. OpenBUGS has to be closed after it finishes running its
# algorithm so the results can return to R.
gumbel.p.min = bugs(data = data.p.min, inits = inits.p.min, 
                    parameters =c("alpha0", "alpha1", "beta", "min_p"),
                    model = "gumbel_p_min.txt",
                    n.thin = 4, n.chains = 1,
                    n.burnin = 5000, n.iter = 12000, codaPkg=F)

# Convergence tests
gumbel.p.min.mcmc <- as.mcmc(gumbel.p.min$sims.matrix[, ])
geweke.diag(gumbel.p.min.mcmc)
raftery.diag(gumbel.p.min.mcmc)
heidel.diag(gumbel.p.min.mcmc)

# Results table, with parameters and weight estimates
hpd.interval.p.min <- HPDinterval(gumbel.p.min.mcmc)
results.p.min <- data.frame(
  "Estimate"=gumbel.p.min$summary[, 1],
  "HPD lower limit"=hpd.interval.p.min[, 1],
  "HPD upper limit"=hpd.interval.p.min[, 2])
results.p.min <- round(results.p.min, 3)


# PRELIMINAR FIT MAX
# Inputs for maximum
max.p.sample <- prices.max[1:9]
data.p.max <- list(y=max.p.sample, n=length(max.p.sample))
inits.p.max <- function(){list(alpha=rnorm(1), beta=1)}

# Creating a text file to run the "bugs" function in OpenBUGS software
sink("gumbel_p_max.txt")
cat("
    model {
    
    # PRIOR DISTRIBUTIONS
    alpha ~ dnorm(0.0, 0.000001)
    beta ~ dnorm(0.0, 0.000001)
    
    # LIKELIHOOD
    for (i in 1:n) {
    y[i] ~ dgumbel(alpha, beta)
    }
    
    # PREDICTIVE
    for(i in 2:3) {
    max_p[i] <- alpha - beta * log(-log(1 - (1 / i)))
    }
    
    }", fill=TRUE)
sink()

# Distribution fit. On Windows, OpenBUGS will open to run the "bugs"
# function. OpenBUGS has to be closed after it finishes running its
# algorithm so the results can return to R.
gumbel.p.max = bugs(data = data.p.max, inits = inits.p.max, 
                    parameters = c("alpha", "beta", "max_p"),
                    model = "gumbel_p_max.txt",
                    n.thin = 4, n.chains = 1,
                    n.burnin = 5000, n.iter = 12000, codaPkg=F)

# Convergence tests
gumbel.p.max.mcmc <- as.mcmc(gumbel.p.max$sims.matrix[, ])
geweke.diag(gumbel.p.max.mcmc)
raftery.diag(gumbel.p.max.mcmc)
heidel.diag(gumbel.p.max.mcmc)

# Results table, with parameters and weight estimates
hpd.interval.p.max <- HPDinterval(gumbel.p.max.mcmc)
results.p.max <- data.frame(
  "Estimate"=gumbel.p.max$summary[, 1],
  "HPD lower limit"=hpd.interval.p.max[, 1],
  "HPD upper limit"=hpd.interval.p.max[, 2])
results.p.max <- round(results.p.max, 3)


# UPDATED FIT MIN
# Inputs for minimum
trend <- seq(1, 12)
min.sample <- prices.min
data.min <- list(x=trend, y=-min.sample, n=length(min.sample))
inits.min <- function(){ list(alpha0=-50, beta=2, alpha1=0.1)}

# Creating a text file to run the "bugs" function in OpenBUGS software
sink("gumbel_min.txt")
cat("
    model {
    
    # PRIOR DISTRIBUTIONS
    alpha0 ~ dnorm(0.0, 0.000001)
    beta ~ dnorm(0.0, 0.000001)
    alpha1 ~ dnorm(0.0, 0.0001)
    
    # LIKELIHOOD
    for (i in 1:n) {
    y[i] ~ dgumbel(alpha[i], beta)
    alpha[i] <- alpha0 + alpha1 * x[i]
    }
    
    # PREDICTIVE
    for(i in 2:6) {
    min_p[i] <- ((alpha0 + alpha1 * (12 + i)) -
                 beta * log(-log(1 - (1 / i))))
    }
    
    }", fill=TRUE)
sink()

# Distribution fit. On Windows, OpenBUGS will open to run the "bugs"
# function. OpenBUGS has to be closed after it finishes running its
# algorithm so the results can return to R.
gumbel.min = bugs(data = data.min, inits = inits.min, 
                  parameters =c("alpha0", "alpha1", "beta", "min_p"),
                  model = "gumbel_min.txt",
                  n.thin = 4, n.chains = 1,
                  n.burnin = 5000, n.iter = 12000, codaPkg=F)

# Convergence tests
gumbel.min.mcmc <- as.mcmc(gumbel.min$sims.matrix[, ])
geweke.diag(gumbel.min.mcmc)
raftery.diag(gumbel.min.mcmc)
heidel.diag(gumbel.min.mcmc)

# Results table, with parameters and weight estimates
hpd.interval.min <- HPDinterval(gumbel.min.mcmc)
results.min <- data.frame(
  "Estimate"=gumbel.min$summary[, 1],
  "HPD lower limit"=hpd.interval.min[, 1],
  "HPD upper limit"=hpd.interval.min[, 2])
results.min <- round(results.min, 3)

# UPDATED FIT MAX 
# Inputs for maximum
max.sample <- prices.max
data.max <- list(y=max.sample, n=length(max.sample))
inits.max <- function(){list(alpha=rnorm(1), beta=1)}

# Creating a text file to run the "bugs" function in OpenBUGS software
sink("gumbel_max.txt")
cat("
    model {
    
    # PRIOR DISTRIBUTIONS
    alpha ~ dnorm(0.0, 0.000001)
    beta ~ dnorm(0.0, 0.000001)
    
    # LIKELIHOOD
    for (i in 1:n) {
    y[i] ~ dgumbel(alpha, beta)
    }
    
    # PREDICTIVE
    for(i in 2:6) {
    max_p[i] <- alpha - beta * log(-log(1 - (1 / i)))
    }
    
    }", fill=TRUE)
sink()

# Distribution fit. On Windows, OpenBUGS will open to run the "bugs"
# function. OpenBUGS has to be closed after it finishes running its
# algorithm so the results can return to R.
gumbel.max = bugs(data = data.max, inits = inits.max, 
                  parameters = c("alpha", "beta", "max_p"),
                  model = "gumbel_max.txt",
                  n.thin = 4, n.chains = 1,
                  n.burnin = 5000, n.iter = 12000, codaPkg=F)

# Convergence tests
gumbel.max.mcmc <- as.mcmc(gumbel.max$sims.matrix[, ])
geweke.diag(gumbel.max.mcmc)
raftery.diag(gumbel.max.mcmc)
heidel.diag(gumbel.max.mcmc)

# Results table, with parameters and weight estimates
hpd.interval.max <- HPDinterval(gumbel.max.mcmc)
results.max <- data.frame(
  "Estimate"=gumbel.max$summary[, 1],
  "HPD lower limit"=hpd.interval.max[, 1],
  "HPD upper limit"=hpd.interval.max[, 2])
results.max <- round(results.max, 3)

#### ================== POSTERIOR DENSITY GRAPHS ================== ####
# Data
volume.6 <- schumacher.vol.mcmc[, 9]
price.min.6 <- -gumbel.min.mcmc[, 8]
price.max.6 <- gumbel.max.mcmc[, 7]
price.min.max.6 <- melt(list("min"=price.min.6, "max"=price.max.6))

# Posterior density for volume at 6 years
vol.posterior <- ggplot(
  melt(list("vol"=volume.6)), aes(value, fill=L1)) +
  geom_histogram(
    aes(y=..density..), colour="black", binwidth=1, size=0.3) +
  scale_fill_manual(breaks=c("vol"), values=c("white")) +
  scale_x_continuous(breaks=seq(305, 360, 11), limits=c(305, 360)) +
  scale_y_continuous(breaks=seq(0, 0.5, 0.1), limits=c(0, 0.32)) +
  labs(title="Volume", x=bquote("Volume"~(m^3~ha^-1)), y="p(x)") +
  theme_gray(base_size=12) +
  theme(axis.text=element_text(size=12, color="black"),
        plot.title=element_text(face="bold"),
        legend.position="",
        axis.ticks=element_line(color="black"))

# Posterior density for minimum and maximum price at 6 years
price.posterior <- ggplot(
  price.min.max.6, aes(value, fill=L1)) +
  geom_histogram(
    aes(y=..density..),
    colour="black", binwidth=1, size=0.3, position="stack") +
  scale_fill_manual(
    breaks=c("min", "max"),
    values=c("white", "white")) +
  scale_x_continuous(breaks=seq(12, 72, 10), limits=c(12, 72)) +
  labs(title="Minimum and maximum prices",
       x=bquote("Price"~(R*"$"~m^3)), y=element_blank()) +
  theme_gray(base_size=12) +
  theme(axis.text=element_text(size=12, color="black"),
        plot.title=element_text(face="bold"),
        legend.position="",
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank())

# Joining graphs
posterior.vol.min.max <- plot_grid(vol.posterior, price.posterior)

# Show graphs
posterior.vol.min.max

#### ====================== ECONOMIC ANALYSIS ===================== ####
# NPV calculation
# Data
volume.6 <- schumacher.vol.mcmc[, 9]
price.min.6 <- -gumbel.min.mcmc[, 8]
price.max.6 <- gumbel.max.mcmc[, 7]
disc.rate <- 0.09

# Revenue
rev.min <- sort(
  melt(as.numeric(volume.6) %o% as.numeric(price.min.6))$value)
rev.max <- sort(
  melt(as.numeric(volume.6) %o% as.numeric(price.max.6))$value)
HPDinterval(as.mcmc(rev.min))
HPDinterval(as.mcmc(rev.max))

# Costs
costs <- data.frame(
  age=c(      0,      1,      2,     3,     4,     5,     6),
  val=c(1520.26, 836.82, 546.47, 75.52, 75.52, 75.52, 75.52))

Cj <- c()
for (i in 1 : length(costs$age)) {
  Cj[i] <- costs$val[i] / (1 + disc.rate) ^ costs$age[i]
}

C <- sum(Cj)

# NPV
NPV.min <- (rev.min / (1 + disc.rate) ^ 6) - C
NPV.max <- (rev.max / (1 + disc.rate) ^ 6) - C
HPDinterval(as.mcmc(NPV.min))
HPDinterval(as.mcmc(NPV.max))


# GRAPHS
# Thinning on NPV so the plot function does not overload
index <- c(
  seq(1, length(rev.min), sqrt(length(rev.min))), length(rev.min))

NPV.min.thin <- c()
for (i in 1 : length(index)) {
  NPV.min.thin[i] <- NPV.min[index[i]]
}

NPV.max.thin <- c()
for (i in 1 : length(index)) {
  NPV.max.thin[i] <- NPV.max[index[i]]
}

# Data for graph
min.max.melt <- melt(list("NPV_min"=NPV.min.thin,
                          "NPV_max"=NPV.max.thin))

# Densities for min and max NPV
min.max.dens <- ggplot(
  min.max.melt, aes(value, fill=L1)) +
  geom_histogram(
    aes(y=..count../sum(..count..) * 2),
    colour="black", bins=70, size=0.3, position="stack") +
  scale_fill_manual(
    breaks=c("NPV_min", "NPV_max"),
    values=c("white", "white")) +
  scale_x_continuous(breaks=seq(0, 11000, 2750), limits=c(0, 11000)) +
  labs(title="NPV for minimum and maximum prices",
       x=bquote("NPV"~(R*"$"~ha^-1)), y="p(x)") +
  theme_gray(base_size=12) +
  theme(axis.text=element_text(size=12, color="black"),
        plot.title=element_text(face="bold"),
        legend.position="",
        axis.ticks=element_line(color="black"))

# Show graph
min.max.dens
