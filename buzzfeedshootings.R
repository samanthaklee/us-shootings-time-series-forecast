library(MASS)
library(GeneCycle)
library(TSA)
library(forecast)
library(astsa)
library(tseries)
library(lubridate)
library(dplyr)
library(ggplot2)

setwd("~/Desktop/PSTAT274Project")
#read in the data and set column names for consistency 
data2013 <- read.csv("2013MASTER.csv")
colnames(data2013) <- c("Reported", "Date", "Shooter", "Killed", "Wounded", "Location", "Article1",
                        "Article2", "Article3", "Article4", "Article5", "Article6", "Article7", "Article8",
                        "Article9")
data2014 <- read.csv("2014MASTER.csv")
colnames(data2014) <- c("Reported", "Date", "Shooter", "Killed", "Injured", "Location", "Article", "Article1",
                        "Article2", "Article3", "Article4", "Article5")
data2015 <- read.csv("CURRENT2015.csv")
data2015$X <- NULL #don't want this
colnames(data2015) <- c("Date", "Shooter", "Killed", "Injured", "Location", "Article", "Article1", "Article2",
                        "Article3", "Article4")

#function to merge all 3 dataframes in R 
MyMerge <- function(x, y){
  df <- merge(x, y, all.x= TRUE, all.y= TRUE)
  return(df)
}
dataTotal <- Reduce(MyMerge, list(data2013, data2014, data2015))

#arrange the entire dataset by date 
dataTotal <- dataTotal %>% mutate(Date = as.Date(Date, "%m/%d/%Y")) %>%
  arrange(Date)

dataTotal$Reported[is.na(dataTotal$Reported)] <- 1 #create a count column for each data point
dataTotal$Date <- as.Date(dataTotal$Date, format="%d/%m/%Y")

# Aggregate over week number and number killed
dataTotal$weeknos <- (interval(min(dataTotal$Date), dataTotal$Date) %/% weeks(1)) + 1
buzzfeed_killed <- aggregate(Killed~weeknos, FUN=sum, data=dataTotal, na.rm=TRUE)
buzzfeed_killed <- data.frame(buzzfeed_killed)

#Remove last 12 rows for forecasting purposes 
buzzfeed_df_removed <- buzzfeed_killed[-c(140:152), ]
View(buzzfeed_df_removed)

#time series
buzzfeed_ts <- ts(buzzfeed_df_removed$Killed) #nonstationary!!
plot(buzzfeed_ts, ylab = "2013 - 2015 Shootings in the US", main = "Weekly Killed in 2013 - 2015 US Shootings")
#linMod <- lm(buzzfeed_ts~as.numeric(1:length(buzzfeed_ts)))
#abline(linMod, col = "red")

#variance to reference for volatility
var(buzzfeed_ts)

#The Augmented Dickey–Fuller (ADF) t-statistic test: 
#small p-values suggest the data is stationary and 
#doesn’t need to be differenced stationarity.
adf.test(buzzfeed_ts, alternative = "stationary")

#The Ljung-Box test examines whether there is significant evidence for non-zero 
#correlations at lags 1-20. Small p-values (i.e., less than 0.05) 
#suggest that the series is stationary.
Box.test(buzzfeed_ts, type = "Ljung-Box")

#The Kwiatkowski-Phillips-Schmidt-Shin (KPSS) test; 
#here accepting the null hypothesis means that the series is stationarity, 
#and small p-values suggest that the series is not stationary and a 
#differencing is required.
kpss.test(buzzfeed_ts)

########seasonal decomposition to find the additive trend, seasonal, irregular components 
decomp <- stl(buzzfeed_ts, s.window = "period")
plot(decomp)

#Take difference at lag 1 twice. Variance of diff2 is greater. 
#Overdifferencing according to the variances
buzzfeed_diff1 <- diff(buzzfeed_ts, lag = 1)
plot(buzzfeed_diff1, ylab = "US Shootings in 2013 - 2015", main = "First Difference on Lag 1")
#linMod2 <- lm(buzzfeed_diff1~as.numeric(1:length(buzzfeed_diff1)))
#abline(linMod, col = "red")

#second difference at lag 1
buzzfeed_diff2 <- diff(buzzfeed_diff1, lag = 1)
plot(buzzfeed_diff2, ylab = "US Shootings in 2013 - 2015", main = "Second Difference on Lag 1")

#compare the variances 
var(buzzfeed_diff1)
var(buzzfeed_diff2)

#Run stationarity tests
Box.test(buzzfeed_diff1, type = "Ljung-Box")
adf.test(buzzfeed_diff1, alternative = "stationary")

#decomposition
diff_decomp <- stl(buzzfeed_diff1, s.window = "periodic")
plot(diff_decomp)

#ACF and PACF of differenced and transformed data
par(mfrow=c(1,2))
acf(buzzfeed_diff1)
acf(buzzfeed_diff1, type = "partial")

#automated arima forecast
auto_arima <- auto.arima(buzzfeed_diff1, trace = TRUE)
summary(auto_arima)
BIC(auto_arima)

#try out our own models
model1 <- arima(buzzfeed_diff1, order = c(4,0,1))
summary(model1)
BIC(model1)

model2 <- arima(buzzfeed_diff1, order = c(1,0,14))
summary(model2)
BIC(model2)

model3 <- arima(buzzfeed_diff1, order = c(1,0,2))
summary(model3)
BIC(model3)

model4 <- arima(buzzfeed_diff1, order = c(2,0,14))
summary(model4)
BIC(model4)

model5 <- arima(buzzfeed_diff1, order = c(4,1,2))
summary(model5)
BIC(model5)

model6 <- arima(buzzfeed_diff1, order = c(4,1,14))
summary(model6)
BIC(model6)

model7 <- arima(buzzfeed_diff1, order = c(3,1,2))
summary(model7)
BIC(model7)

model8 <- arima(buzzfeed_diff1, order = c(1,1,1))
summary(model8)
BIC(model8)

model9 <- arima(buzzfeed_diff1, order = c(2,1,4))
summary(model9)
BIC(model9)

model10 <- arima(buzzfeed_diff1, order = c(4,1,2))
summary(model10)
BIC(model10)

#we will use model 3, model 1, auto.arima, model 9

#MODEL A DIAGNOSTICS
plot(model3$residuals, main="Model A Residuals", ylab = "Residuals")
shapiro.test(model3$residuals)

qqnorm(model3$residuals)

paste("Box-Pierce")
Box.test(model3$residuals, lag = 5,  type = "Box-Pierce", fitdf=2)
paste("Ljung-Box") 
Box.test(model3$residuals, lag = 5, type = "Ljung-Box", fitdf=2)
paste("McLeod-Li")
Box.test((model3$residuals)^2, lag=5, type="Ljung-Box")

par(mfrow=c(1,1))
acf(model7$residuals, na.action=na.pass, main = "ACF of Residuals for Model A")
pacf(model7$residuals, na.action = na.pass, main = "PACF of Residuals for Model A")

#MODEL B DIAGNOSTICS
plot(model1$residuals, main="Model B Residuals", ylab = "Residuals")
shapiro.test(model1$residuals)
qqnorm(model1$residuals)

acf(model1$residuals, na.action=na.pass, main = "ACF of Residuals for Model B")
pacf(model1$residuals, na.action=na.pass, main = "PACF of Residuals for Model B")

paste("Box-Pierce")
Box.test(model1$residuals, lag = 5, type = "Box-Pierce", fitdf=2)
paste("Ljung-Box")
Box.test(model1$residuals, lag = 5, type = "Ljung-Box", fitdf=2)
paste("McLeod-Li")
Box.test((model1$residuals)^2, lag=5, type="Ljung-Box")

#MODEL C DIAGNOSTICS
plot(auto_arima$residuals, main="Model C Residuals", ylab = "Residuals")
shapiro.test(auto_arima$residuals)
qqnorm(auto_arima$residuals)
acf(auto_arima$residuals, na.action=na.pass, main = "ACF of Residuals for Model C")
pacf(auto_arima$residuals, na.action=na.pass, main = "PACF of Residuals for Model C")

paste("Box-Pierce")
Box.test(auto_arima$residuals, lag = 5, type = "Box-Pierce", fitdf=2)
paste("Ljung-Box")
Box.test(auto_arima$residuals, lag = 5, type = "Ljung-Box", fitdf=2)
paste("McLeod-Li")
Box.test((auto_arima$residuals)^2, lag=5, type="Ljung-Box")


#MODEL D DIAGNOSTICS
plot(model9$residuals, main="Model D Residuals", ylab = "Residuals")
shapiro.test(model9$residuals)
qqnorm(model9$residuals)
acf(model9$residuals, na.action=na.pass, main = "ACF of Residuals for Model D")
pacf(model9$residuals, na.action=na.pass, main = "PACF of Residuals for Model D")

paste("Box-Pierce")
Box.test(model9$residuals, lag = 5, type = "Box-Pierce", fitdf=2)
paste("Ljung-Box")
Box.test(model9$residuals, lag = 5, type = "Ljung-Box", fitdf=2)
paste("McLeod-Li")
Box.test((model9$residuals)^2, lag=5, type="Ljung-Box")



#Begin forecasting 
#match forecast to our removed data
#the last twelve points of our data to compare against the forecast
buzzfeed_last_12 <- buzzfeed_killed[-c(1:139),]
View(buzzfeed_last_12)
orig_removed <- buzzfeed_last_12$Killed
View(orig_removed)

#fitted model C
fit <- arima(buzzfeed_ts, order= c(1,0,1), method = "ML", xreg = 1:length(buzzfeed_ts))
buzzfeed_predict <- predict(fit, n.ahead = 13, newxreg = (length(buzzfeed_ts)+1): length(buzzfeed_ts)+13)

#confidence interval bounds
upper <- buzzfeed_predict$pred + 2*buzzfeed_predict$se
lower <- buzzfeed_predict$pred - 2*buzzfeed_predict$se

#plot forecast
plot(buzzfeed_ts, xlim = c(100, 160), ylim = c(min(lower), max(upper)), ylab = "2013 to 2015 Shootings in the US",
     main = "Forecast of Shootings in the US")
lines(upper, col = "blue", lty = "dashed")
lines(lower, col = "blue", lty = "dashed")
points(140:152, buzzfeed_predict$pred, col = "green")
points(140:152,orig_removed, pch = "*")

#fitted model D
fit2 <- arima(buzzfeed_ts, order= c(2,1,4), method = "ML", xreg = 1:length(buzzfeed_ts))
buzzfeed_predict2 <- predict(fit2, n.ahead = 13, newxreg = (length(buzzfeed_ts)+1): length(buzzfeed_ts)+13)

#confidence interval bounds
upper2 <- buzzfeed_predict2$pred + 2*buzzfeed_predict2$se
lower2 <- buzzfeed_predict2$pred - 2*buzzfeed_predict2$se

#plot forecast
plot(buzzfeed_ts, xlim = c(100, 160), ylim = c(min(lower2), max(upper2)), ylab = "2013 to 2015 Shootings in the US",
     main = "Forecast of Shootings in the US")
lines(upper2, col = "blue", lty = "dashed")
lines(lower2, col = "blue", lty = "dashed")
points(140:152, buzzfeed_predict2$pred, pch = "*")
points(140:152,orig_removed, col = "green")

#coefficients
summary(fit)
summary(fit2)


#spectral analysis
periodogram(buzzfeed_ts) 

#smoothed periodogram
del<-0.1 # sampling interval
x.spec <- spectrum(buzzfeed_ts,log="no",span=10,plot=FALSE)
spx <- x.spec$freq/del #cycles per unit of time
spy <- 2*x.spec$spec # multiply spectral density by 2 to match variance of time series
plot(spy~spx,xlab="Frequency",ylab="Spectral density",type="l", main = "Smoothed Periodogram")

#fisher test
fisher.g.test(residuals(fit))

#ks test
cpgram(residuals(fit), main = "")


#buzzfeed_period <-spec.pgram(buzzfeed_ts,fast=FALSE, taper =0.0);


############################################################################################################
#unused portion in the report

set.seed(42)
bcTransform <- boxcox(buzzfeed_ts ~ as.numeric(1:length(buzzfeed_ts)))
trans <- bcTransform$x[which.max(bcTransform$y)]
trans
#log is not the best transformation 

buzzLog <- log(buzzfeed_ts)
plot.ts(buzzLog)

class(buzzLog) #ts
#differencing of the log transformation
logts_shooting <- diff(buzzLog)
plot(logts_shooting)
acf(logts_shooting)
acf(logts_shooting, type = "partial")
var(logts_shooting) #0.088

#try differencing on original data
#buzzfeed_ds_diff_OLD <- diff(buzzfeed_ts)
#plot(buzzfeed_ds_diff_OLD, ylab = "Differenced Data of US Shootings", main = "Differenced Data Plot")

#run stationary tests
#Box.test(buzzfeed_ds_diff_OLD, type = "Ljung-Box")
#adf.test(buzzfeed_ds_diff_OLD, alternative = "stationary")

#difference square root data
#buzzfeed_ds_diff <- diff(sqrt_trans)
#buzzfeed_ds_diff <- ts(buzzfeed_ds_diff, start = c(2013, 1), frequency =12)
#plot(buzzfeed_ds_diff, ylab = "Square Root Difference of Data", main = "Square Root Differenced Data of US Shootings")

#stationary tests
#Box.test(buzzfeed_ds_diff, type = "Ljung-Box")
#adf.test(buzzfeed_ds_diff, alternative = "stationary")
#var(buzzfeed_ds_diff)

#forecast function
model7_forecast <- forecast(model7)
plot(model7_forecast)
accuracy(model7_forecast)

#spectral analysis
par(mfcol = c(2,2))
plot.ts(ts_shootings, main = "Differenced Time Series")
mvspec(ts_shootings, spans = c(5,5), plot = TRUE) #nonparametric spectral estimate
arma.spec(ar = -0.7608, ma = 0.7152)


freq3 = 0.496402878
A3 = rnorm(1, 0, 4)
B3 = rnorm(1,0,4)
z = 2*pi*buzzfeed_ts
x = A3*cos(z*freq3) + B3*sin(z*freq3)  + rnorm(96,0,1)
plot(buzzfeed_ts,x,type='o',ylab=expression(x[t]))
periodogram(buzzfeed_ts); abline(h=0)



freq1=0.100719424
freq2=0.424460432 # freq1=0.46875; freq2=0.04166667
A1=rnorm(1,0,2) 
B1=rnorm(1,0,2) # random amplitudes
A2=rnorm(1,0,3)
B2=rnorm(1,0,3)
w=2*pi*buzzfeed_ts
y=A1*cos(w*freq1) + B1*sin(w*freq1) + A2*cos(w*freq2) + B2*sin(w*freq2) + rnorm(96,0,1)
plot(buzzfeed_ts,y,type='o',ylab=expression(y[t]))
periodogram(buzzfeed_ts); abline(h=0)

#fitted model A
fit3 <- arima(buzzfeed_ts, order= c(3,1,2), method = "ML", xreg = 1:length(buzzfeed_ts))
buzzfeed_predict3 <- predict(fit3, n.ahead = 13, newxreg = (length(buzzfeed_ts)+1): length(buzzfeed_ts)+13)

#confidence interval bounds
upper3 <- buzzfeed_predict3$pred + 2*buzzfeed_predict3$se
lower3 <- buzzfeed_predict3$pred - 2*buzzfeed_predict3$se

#plot forecast
plot(buzzfeed_ts, xlim = c(100, 160), ylim = c(min(lower3), max(upper3)), ylab = "2013 to 2015 Shootings in the US",
     main = "Forecast of Shootings in the US")
lines(upper3, col = "blue", lty = "dashed")
lines(lower3, col = "blue", lty = "dashed")
points(140:152, buzzfeed_predict3$pred, pch = "*")
points(140:152,orig_removed, col = "green")


#fitted model 10
fit4 <- arima(buzzfeed_ts, order= c(4,1,2), method = "ML", xreg = 1:length(buzzfeed_ts))
buzzfeed_predict4 <- predict(fit4, n.ahead = 13, newxreg = (length(buzzfeed_ts)+1): length(buzzfeed_ts)+13)

#confidence interval bounds
upper4 <- buzzfeed_predict4$pred + 2*buzzfeed_predict4$se
lower4 <- buzzfeed_predict4$pred - 2*buzzfeed_predict4$se

#plot forecast
plot(buzzfeed_ts, xlim = c(100, 160), ylim = c(min(lower3), max(upper3)), ylab = "2013 to 2015 Shootings in the US",
     main = "Forecast of Shootings in the US")
lines(upper4, col = "blue", lty = "dashed")
lines(lower4, col = "blue", lty = "dashed")
points(140:152, buzzfeed_predict4$pred, pch = "*")
points(140:152,orig_removed, col = "green")
