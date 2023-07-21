library(ggplot2)
library(ADMUR)
library(DEoptimR)

load(file = here::here("analysis/data/derived_data/shore_pd_models.RData"))

SPD <- as.data.frame( rowSums(pd) )
# normalise
SPD <- SPD/( sum(SPD) * 5 )
plotPD(SPD)

library(pso)

normalize <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}

mmpd <- as.data.frame(apply(pd, 2, normalize))

CPL.1pso <- psoptim(par = NA,
                    lower = 0,
                    upper = 0.999999,
                    fn = objectiveFunction,
                    PDarray = pd, type = 'CPL',
                    control = list(trace = 1))
CPL.2pso <- psoptim(par = rep(NA, 3),
                    lower = rep(0, 3),
                    upper = rep(0.999999, 3),
                    fn = objectiveFunction,
                    PDarray = pd,
                    type = 'CPL',
                    control = list(trace = 1))
CPL.3pso <- psoptim(par = rep(NA, 5),
                    lower = rep(0, 5),
                    upper = rep(0.999999, 5),
                    fn = objectiveFunction,
                    PDarray = pd,
                    type = 'CPL',
                    control = list(trace = 1))

CPL1pso <- convertPars(pars=CPL.1pso$par, years=minage:maxage, type='CPL')
CPL2pso <- convertPars(pars=CPL.2pso$par, years=minage:maxage, type='CPL')
CPL3pso <- convertPars(pars=CPL.3pso$par, years=minage:maxage, type='CPL')

ggplot() +
  geom_bar(aes(x = as.numeric(rownames(SPD)),
               SPD[,1]),
           stat = "identity", col = "grey") +
  geom_line(aes(x = as.numeric(rownames(SPD)),
                y = SPD[,1])) +
  geom_line(data = CPL1pso, aes(year, pdf),
            linewidth = 1) +
  geom_line(data = CPL2pso, aes(year, pdf), col = "red",
            linewidth = 1) +
  geom_line(data = CPL3pso, aes(year, pdf), col = "blue",
            linewidth = 1)




round(ncol(pd) / 6)

s1 <- pd[, sample(ncol(pd), 6)]
s1spd <- as.data.frame( rowSums(s1) )
s1spd <- s1spd/( sum(s1) * 5 )
plotPD(s1spd)

unifs <- -objectiveFunction(pars = NULL, PDarray = s1, type = 'uniform')
cpl_1s <- JDEoptim(lower = 0, upper = 1, fn = objectiveFunction, PDarray = s1,
                  type = 'CPL', NP = 20, trace = TRUE)
cpl_2s <- JDEoptim(lower = rep(0, 3), upper = rep(1, 3), fn = objectiveFunction,
                  PDarray = s1, type = 'CPL', NP = 60,  maxiter = 400 * 3,
                  trace = TRUE)
cpl_3s <- JDEoptim(lower = rep(0, 5), upper = rep(1, 5), fn = objectiveFunction,
                  PDarray = s1, type = 'CPL', NP = 100, maxiter = 400 * 10,
                  trace = TRUE)
cpl_4s <- JDEoptim(lower = rep(0, 7), upper = rep(1, 7), fn = objectiveFunction,
                  PDarray = s1, type = 'CPL', NP = 140, maxiter = 400 * 10,
                  trace = TRUE)
cpl_5s <- JDEoptim(lower = rep(0, 9), upper = rep(1, 9), fn = objectiveFunction,
                  PDarray = s1, type='CPL', NP = 180, maxiter = 400 * 20,
                  trace = TRUE)
cpl_6s <- JDEoptim(lower = rep(0, 11), upper = rep(1,11), fn = objectiveFunction,
                  PDarray = s1, type = 'CPL', NP = 220, maxiter =  400 * 30,
                  trace = TRUE)

uni <- convertPars(pars = NULL, years = minage:maxage, type = "uniform")
uni$name <- "Uniform"
cpl1s1 <- convertPars(pars = cpl_1s$par, years = minage:maxage, type = 'CPL')
cpl1s1$name <- "1-CPL"
cpl2s1 <- convertPars(pars = cpl_2s$par, years = minage:maxage, type = 'CPL')
cpl2s1$name <- "2-CPL"
cpl3s1 <- convertPars(pars = cpl_3s$par, years = minage:maxage, type = 'CPL')
cpl3s1$name <- "3-CPL"
cpl4s1 <- convertPars(pars = cpl_4s$par, years = minage:maxage, type = 'CPL')
cpl4s1$name <- "4-CPL"
cpl5s1 <- convertPars(pars = cpl_5s$par, years = minage:maxage, type = 'CPL')
cpl5s1$name <- "5-CPL"
cpl6s1 <- convertPars(pars = cpl_6s$par, years = minage:maxage, type = 'CPL')
cpl6s1$name <- "6-CPL"

s1models <- rbind(uni, cpl1s1, cpl2s1, cpl3s1, cpl4s1, cpl5s1, cpl6s1)

ggplot() +
  geom_line(aes(x = as.numeric(rownames(s1spd)), y = s1spd[,1])) +
  geom_line(data = s1models, aes(x = year, y = pdf, col = name)) +
  scale_x_reverse() +
  labs(x = "BP", y = "PD")

s2 <- pd[, sample(ncol(pd), 32)]
s2spd <- as.data.frame( rowSums(s2) )
s2spd <- s2spd/( sum(s2) * 5 )
plotPD(s2spd)

unifs <- -objectiveFunction(pars = NULL, PDarray = s2, type = 'uniform')
cpl_1s2 <- JDEoptim(lower = 0, upper = 1, fn = objectiveFunction, PDarray = s2,
                   type = 'CPL', NP = 20, trace = TRUE)
cpl_2s2 <- JDEoptim(lower = rep(0, 3), upper = rep(1, 3), fn = objectiveFunction,
                   PDarray = s2, type = 'CPL', NP = 60,  maxiter = 400 * 3,
                   trace = TRUE)
cpl_3s2 <- JDEoptim(lower = rep(0, 5), upper = rep(1, 5), fn = objectiveFunction,
                   PDarray = s2, type = 'CPL', NP = 100, maxiter = 400 * 10,
                   trace = TRUE)
cpl_4s2 <- JDEoptim(lower = rep(0, 7), upper = rep(1, 7), fn = objectiveFunction,
                   PDarray = s2, type = 'CPL', NP = 140, maxiter = 400 * 10,
                   trace = TRUE)
cpl_5s2 <- JDEoptim(lower = rep(0, 9), upper = rep(1, 9), fn = objectiveFunction,
                   PDarray = s2, type='CPL', NP = 180, maxiter = 400 * 20,
                   trace = TRUE)
cpl_6s2 <- JDEoptim(lower = rep(0, 11), upper = rep(1,11), fn = objectiveFunction,
                   PDarray = s2, type = 'CPL', NP = 220, maxiter =  400 * 30,
                   trace = TRUE)

uni <- convertPars(pars = NULL, years = minage:maxage, type = "uniform")
uni$name <- "Uniform"
cpl1s2 <- convertPars(pars = cpl_1s2$par, years = minage:maxage, type = 'CPL')
cpl1s2$name <- "1-CPL"
cpl2s2 <- convertPars(pars = cpl_2s2$par, years = minage:maxage, type = 'CPL')
cpl2s2$name <- "2-CPL"
cpl3s2 <- convertPars(pars = cpl_3s2$par, years = minage:maxage, type = 'CPL')
cpl3s2$name <- "3-CPL"
cpl4s2 <- convertPars(pars = cpl_4s2$par, years = minage:maxage, type = 'CPL')
cpl4s2$name <- "4-CPL"
cpl5s2 <- convertPars(pars = cpl_5s2$par, years = minage:maxage, type = 'CPL')
cpl5s2$name <- "5-CPL"
cpl6s2 <- convertPars(pars = cpl_6s2$par, years = minage:maxage, type = 'CPL')
cpl6s2$name <- "6-CPL"

s2models <- rbind(uni, cpl1s2, cpl2s2, cpl3s2, cpl4s2, cpl5s2, cpl6s2)

ggplot() +
  geom_line(aes(x = as.numeric(rownames(s2spd)), y = s2spd[,1])) +
  geom_line(data = s2models, aes(x = year, y = pdf, col = name)) +
  scale_x_reverse() +
  labs(x = "BP", y = "PD")


sampsize <- c(6, 32, 60, 120, 240, 480)


for(i in 1:sampsize){
  tmppd <- pd[, sample(ncol(pd), sampsize[i])]

  cpl_4tmp <- JDEoptim(lower = rep(0, 7),
                       upper = rep(1, 7),
                       fn = objectiveFunction,
                      PDarray = tmppd,
                      type = 'CPL',
                      NP = 140)

  cpl_5tmp <- JDEoptim(lower = rep(0, 9),
                       upper = rep(1, 9),
                       fn = objectiveFunction,
                      PDarray = tmppd,
                      type='CPL',
                      NP = 180,
                      maxiter = 400 * 9,
                      trace = TRUE)

  cpl_6tmp <- JDEoptim(lower = rep(0, 11),
                       upper = rep(1,11),
                       fn = objectiveFunction,
                       PDarray = tmppd,
                       type = 'CPL',
                       NP = 220,
                       maxiter =  400 * 11,
                       trace = TRUE)

  print(paste("Sample size =", sampsize[i]))
}

objectiveFunction <- function(pars, PDarray, type, taphonomy=FALSE){

  if(!is.data.frame(PDarray))stop('PDarray must be a data frame')
  years <- as.numeric(row.names(PDarray))
  model <- convertPars(pars,years,type,taphonomy)
  loglik <- loglik(PDarray, model)

  return(-loglik)}

convertPars <- function(pars, years, type, taphonomy=FALSE){

  # The model must be returned as a PDF. I.e, the total area must sum to 1.

  # sanity checks
  model.choices <- c('CPL','exp','uniform','norm','sine','cauchy','logistic','power')
  if(!type%in%model.choices)stop(paste('Unknown model type. Choose from:',paste(model.choices,collapse=', ')))
  if('data.frame'%in%class(pars))pars <- as.matrix(pars)
  if('integer'%in%class(years))years <- as.numeric(years)
  if(!'numeric'%in%class(years))stop('years must be a numeric vector')

  if('NULL'%in%class(pars) | 'numeric'%in%class(pars)){
    res <- convertParsInner(pars, years, type, taphonomy)
    return(res)
  }

  if(!'numeric'%in%class(pars)){
    N <- nrow(pars)
    C <- length(years)
    res <- as.data.frame(matrix(,N,C))
    names(res) <- years
    for(n in 1:N)res[n,] <- convertParsInner(pars[n,], years, type, taphonomy)$pdf
  }
  return(res)
}

loglik <- function(PD, model){

  # relative likelihood of a perfectly precise date is the model PDF
  # therefore the relative likelihood of a date with uncertainty is an average of the model PDF, weighted by the date probabilities.
  # Numerically this is the scalar product: sum of (model PDF x date PDF).

  years <- as.numeric(row.names(PD))

  # ensure the date ranges exactly match. If not, interpolate model pdf to match PD.
  check <- identical(years,model$year)
  if(!check) model <- interpolate.model.to.PD(PD, model)

  # ensure model PD is provided as a discretised PDF
  inc <- (years[2]-years[1])
  model$pdf <- model$pdf/(sum(model$pdf)*inc)

  # convert the date PD pdfs to discrete PMFs to perform a weighted average
  PMF <- PD * inc

  # likelihoods weighted by the observational uncertainty
  weighted.PD <- PMF * model$pdf

  # sum all possibilities for each date (a calibrated date's probabilities are OR) to give the relative likelihood for each date.
  liks <- colSums(weighted.PD)

  # calculate the overall log lik given all the dates
  loglik <- sum(log(liks))

  if(is.nan(loglik))loglik <- -Inf
  return(loglik)}

convertParsInner <- function(pars, years, type, taphonomy){

  if(taphonomy){
    p <- length(pars)
    model.pars <- pars[0:(p-2)]
    taph.pars <- pars[(p-1):p]
  }
  if(!taphonomy){
    model.pars <- pars
    taph.pars <- c(0,0)
  }

  if(type=='CPL'){
    tmp <- CPLPDF(years,model.pars)
  }

  # incorporate taphonomy
  taph <- (years + taph.pars[1])^taph.pars[2]
  tmp <- tmp * taph
  inc <- years[2]-years[1]
  pdf <- tmp/(sum(tmp)*inc)

  res <- data.frame(year = years, pdf = pdf)
  return(res)}

CPLPDF <- function(x,pars){
  hinges <- CPLparsToHinges(pars, x)
  if(any(is.na(hinges$pdf))){
    pdf <- NA
  } else {
    pdf <- approx(x=hinges$year, y=hinges$pdf, xout=x)$y
  }
  return(pdf)}

CPLparsToHinges <- function(pars, years){

  if('numeric'%in%class(pars)){
    res <- CPLparsToHingesInner(pars, years)
    return(res)}

  N <- nrow(pars)
  C <- (ncol(pars)+1)/2 +1
  yr <- pdf <- as.data.frame(matrix(,N,C))
  names(yr) <- paste('yr',1:C,sep='')
  names(pdf) <- paste('pdf',1:C,sep='')
  for(n in 1:N){
    x <- CPLparsToHingesInner(pars[n,],years)
    yr[n,] <- x$year
    pdf[n,] <- x$pdf
  }
  res <- cbind(yr,pdf)
  return(res)}

CPLparsToHingesInner <- function(pars, years){

  # must be odd, as (2n-1 parameters where n=number of pieces)
  cond <- ((length(pars)+1) %% 2) == 0
  if(!cond)stop('A CPL model must have an odd number of parameters')

  # parameters must be between 0 and 1
  if(sum(pars>1 | pars<0)!=0)stop('CPL parameters must be between 0 and 1')

  if(length(pars)==1){
    x.par <- c()
    y.par <- pars
  }

  if(length(pars)!=1){
    x.par <- pars[1:((length(pars)-1)/2)]
    y.par <- pars[(length(x.par)+1):length(pars)]
  }

  # conversion of pars to raw hinge coordinates x between 0 and 1,  and y between 0 and Inf
  # much more efficient stick breaking algorithm for x
  # mapping for y (0 to 1) -> (0 to Inf) using (1/(1-y)^2)-1
  # y0 is arbitrarily fixed at 3 since (1/(1-0.5)^2)-1
  xn <- length(x.par)
  if(xn>=1)proportion <- qbeta(x.par, 1 , xn:1)
  if(xn==0)proportion <- c()
  x.raw <- c(0,1-cumprod(1 - proportion),1)
  y.raw <- c(3, (1/(1-y.par)^2)-1)

  # convert x.raw from 0 to 1, to years
  x <- x.raw * (max(years)-min(years)) + min(years)

  # area under curve
  widths <- diff(x)
  mids <- 0.5*(y.raw[1:(xn+1)]+y.raw[2:(xn+2)])
  area <- sum(widths*mids)

  # convert y.raw to PD
  y <- y.raw/area

  # store
  d <- data.frame(year=x, pdf=y)
  return(d)}
