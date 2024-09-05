library(xts) #Generate Time Series
library(PerformanceAnalytics) #Return calculation
library(MASS) #To get the inverse of a matrix
library(timeDate) #To identify Weekends
library(quantmod) #Get data from yahoo finance


#The objective behind this simulations is forecast the prices of some GROUP stocks and build a portfolio using HRP to estimate future performance. 


#Get data from Yahoo finance

rv <- c('EWC','TSLA','MCHI','AAPL') #Names of stocks to simulate forecast.
startdate <- c('2023-01-01') #Start date of data
enddate <- Sys.Date()
for(i in 1:length(rv)){
  if(i == 1){
    datos <- getSymbols(rv[i], src='yahoo', from = startdate, to = enddate,auto.assign = FALSE)[,4]
  }else{
    datos <- cbind(datos, getSymbols(rv[i], src='yahoo', from = startdate, to = enddate,auto.assign = FALSE)[,4]  )
  }
}
colnames(datos) <- rv
nombres <- rv

#Stage 1: Transform and process data-------------------------------------------------------------------
datos2 <- na.omit(Return.calculate(datos))                      #1.1: Time series returns
medias <- apply(datos2,2,mean)                                  #1.2: Mean vector
varianza <- cov(datos2)                                         #1.2: Covariance matrix
cholesky <- t(chol(varianza))                                   #1.2: Cholesky(L) decomposition
#---------------------------------------------------------------------------------------------------------------


#Sateg 2: Forecast------------------------------------------------------------------------------------------
years <- 1                                                      #2.1: Years to simulate
simulaciones <- 100                                            #2.1: Number of simulations
d <- datos[nrow(datos),]                                        #2.1: Last day prices
#Dates to forecast
fechas <- matrix(index(d),nrow=1,ncol=1)
for(i in 1:(252*years)){
  fecha <- as.Date( fechas[nrow(fechas),] + 1)
  while(isWeekend(fecha) == TRUE){
    fecha <- as.Date(as.numeric(fecha) + 1)
  }
  fechas <- rbind(fechas,fecha)}
fechas <- fechas[-1,]
forecasts <- list()
for(k in 1:simulaciones){
  forecast <- matrix(as.numeric(d),nrow=1,ncol=ncol(datos2))
  for(i in 1:(252*years)){
    randomnumbers <- as.matrix(c(rnorm(ncol(datos2),0,.5)),ncol=1,nrow= ncol(datos2) )  #2.2: Random numbers generation (x)
    v <- t(cholesky%*%randomnumbers) + medias                   #2.3: u + Lx
    forecast <- rbind(forecast,forecast[nrow(forecast),]*(v+1)) #2.4: Transform returns to prices
  }
  forecast <- forecast[-1,]
  forecasts[[k]] <- as.xts(forecast, as.Date(fechas))
}
#---------------------------------------------------------------------------------------------------------------
#Plots-----
plots <- list()
for(v in 1:ncol(datos)){
  instrumento1 <- rbind(datos[,v],forecasts[[1]][,v])
  for(i in 2:simulaciones){
    instrumento1 <- merge(instrumento1,rbind(datos[,v],forecasts[[i]][,v]))
  }
  plots[[v]] <- plot(instrumento1,main=nombres[v])
}
#-----
#Equally weighted portfolio  ----
w <- rep(1/length(rv),length(rv))
wew <- w
capitalinicial <- 1000
for(i in 1:length(forecasts)){
  r <- na.omit(CalculateReturns(forecasts[[i]]))
  for(j in 1:ncol(r)){
    r[,j] <- r[,j]*w[j]
  }
  rportewdiario <- apply(r,1,sum)
  ac <- 1
  for(k in 1:length(rportewdiario)){
    ac <- c(ac, ac[k]*(1+ as.numeric(rportewdiario[k]) ))
  }
  ac <- as.xts(ac, order.by =  index(forecasts[[i]])     )*capitalinicial
  if(i == 1){
    PortafoliosEW <- ac
  }else{
    PortafoliosEW <- cbind(PortafoliosEW,ac)
  }
}
colnames(PortafoliosEW) <- c(1:ncol(PortafoliosEW))

#----

#HRP Portfolio----
getIVP <- function(covMat) {
  invDiag <- 1/diag(as.matrix(covMat))
  weights <- invDiag/sum(invDiag)
  return(weights)
}
getClusterVar <- function(covMat, cItems) {
  covMatSlice <- covMat[cItems, cItems]
  weights <- getIVP(covMatSlice)
  cVar <- t(weights) %*% as.matrix(covMatSlice) %*% weights
  return(cVar)
}
getRecBipart <- function(covMat, sortIx) {
  w <- rep(1,ncol(covMat))
  w <- recurFun(w, covMat, sortIx)
  return(w)
}
recurFun <- function(w, covMat, sortIx) {
  subIdx <- 1:trunc(length(sortIx)/2)
  cItems0 <- sortIx[subIdx]
  cItems1 <- sortIx[-subIdx]
  cVar0 <- getClusterVar(covMat, cItems0)
  cVar1 <- getClusterVar(covMat, cItems1)
  alpha <- 1 - cVar0/(cVar0 + cVar1)
  w[cItems0] <- w[cItems0] * as.numeric(alpha)
  w[cItems1] <- w[cItems1] * (1-as.numeric(alpha))
  if(length(cItems0) > 1) {
    w <- recurFun(w, covMat, cItems0)
  }
  if(length(cItems1) > 1) {
    w <- recurFun(w, covMat, cItems1)
  }
  return(w)
}
#Annualized returns
ps <- c(nrow(datos):253)
retornosanualizados <- c()
for(i in ps){
  retornosanualizados <- rbind(retornosanualizados,  (as.numeric(datos[i])/as.numeric(datos[i-252])) - 1)
}
colnames(retornosanualizados) <- rv
vector_de_retornos_esperados <- apply(retornosanualizados,2,mean)
varianza_de_retornos <- var(retornosanualizados)

#Apply HRP
covMat <- cov(retornosanualizados)
corMat <- cor(retornosanualizados)
clustOrder <- hclust(dist(corMat), method = 'single')$order

w <- getRecBipart(covMat, clustOrder)
names(w) <- rv
whrp <- w
capitalinicial <- 1000
for(i in 1:length(forecasts)){
  r <- na.omit(CalculateReturns(forecasts[[i]]))
  for(j in 1:ncol(r)){
    r[,j] <- r[,j]*w[j]
  }
  rportewdiario <- apply(r,1,sum)
  ac <- 1
  for(k in 1:length(rportewdiario)){
    ac <- c(ac, ac[k]*(1+ as.numeric(rportewdiario[k]) ))
  }
  ac <- as.xts(ac, order.by =  index(forecasts[[i]])     )*capitalinicial
  if(i == 1){
    PortafoliosHRP <- ac
  }else{
    PortafoliosHRP <- cbind(PortafoliosHRP,ac)
  }
}
colnames(PortafoliosHRP) <- c(1:ncol(PortafoliosHRP))

#----


#Compare Portfolios
pew <- data.frame(Portafolio = "EW",
                  Expected_Return = sum(vector_de_retornos_esperados*wew),
                  Variance = t(wew)%*%varianza_de_retornos%*%wew,
                  Diversification = sum(wew**2))
phrp <- data.frame(Portafolio = "HRP",
                   Expected_Return = sum(vector_de_retornos_esperados*whrp),
                   Variance = t(whrp)%*%varianza_de_retornos%*%whrp,
                   Diversification = sum(whrp**2))
portafolios <- rbind(pew,phrp)


#Final Result
plot( PortafoliosEW, main = paste0( "Portafolios EW, Media = ",  round(mean(PortafoliosEW[nrow(PortafoliosEW),]),2)  )  )
plot( PortafoliosHRP, main = paste0( "Portafolios HRP, Media = ",  round(mean(PortafoliosHRP[nrow(PortafoliosHRP),]),2)  )  )
portafolios
