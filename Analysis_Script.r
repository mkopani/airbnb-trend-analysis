#load relevant packages
#Make sure all these packages are installed before running this script
library(ggplot2)
library(forecast)
library(tseries)
library(IDPmisc)
library(zoo)
library(xts)
library(lubridate)
library(dplyr)
library(glmnet)
library(lars)
library(stats)
#####################

#load data
df <- read.csv('data/combined_listing.csv',header=TRUE,na.strings=c("","NA"))

#create a new df to explore price and time data
n <- nrow(df)
price.col <- df[1:n,'price']
price.col = as.numeric(gsub("\\$", "", price.col))
last.review <- as.Date(df[1:n,'last_review'])

#make new dataframe with only price and last review column
ts.df <- data.frame(price.col,last.review)


#order df by date
ts.df <- ts.df[order(as.Date(ts.df$last.review, format="%Y-%m-%d")),]

#better col names and clean up df
colnames(ts.df) <- c('price','date')
ts.df <- NaRV.omit(ts.df)

#get first and last date 
start_date <- ts.df$date[1]
last_date <- ts.df$date[nrow(ts.df)]

#sequence of dates for date column
date <- data.frame(seq(from=start_date,to=last_date,by='days'))
colnames(date) <- c('date')

#remove any data from before 2015. Not needed/too little
ts.df <- subset(ts.df, ts.df$date> '2015-01-01')

#ts plot of raw prices. May be a little too noisy. Some filtering probably needed
plot(as.Date(ts.df$date),ts.df$price, type='l',ylab='price',xlab='Time',col='blue')


#initial plot too noisy. Group by day
meanprice <-  aggregate(ts.df$price,by=list(ts.df$date),mean)
colnames(meanprice) <- c('date','meanprice')
meanprice <- merge(x = date, y = meanprice, by = 'date', all = TRUE)
meanprice$date <- as.Date(meanprice$date)

#remove data from before 2015. Too little of it
meanprice <- na.omit(meanprice[meanprice[['date']] > as.Date('2015-01-01') ,])
meanprice <- data.frame(meanprice)
meanprice$meanprice <-meanprice$meanprice

#plot with aggregated data. ggplot because its prettier.
p3<- ggplot(meanprice,aes(date,meanprice)) + geom_line() +ylab('Daily Average Price') + xlab('Time')
print(p3)


#After consulting Harry Joe, he suggested it might make sense to aggregate data by week and see if there are trends. 
#group by week
daily_xts <- as.xts(meanprice$meanprice,order.by=as.Date(meanprice$date))
weekly_ts <- apply.weekly(daily_xts,mean)


#plot weekly time series
print(plot.xts(weekly_ts,main='Time series of average weekly price'))

#decompose time series
weekly <- ts(as.numeric(weekly_ts), frequency=52)
weekly_components <- decompose(weekly)

#plot decomposistion
plot(weekly_components)

#Holt-Winter exponential smoothing
exposmoothing_forecasts <- HoltWinters(weekly,seasonal = 'additive')
plot(exposmoothing_forecasts)

#testing stationarity
adf.test(weekly,alternative = "stationary")

#Variogram defintion: 
#G(k) = Var(vt+k - vt)/Var(zt+1 - zt)
Var <- var(diff(weekly,lag=1))
Gk <- vector(mode='numeric',length=48)
for (k in 1:48){
  Gk[k] <- var(diff(weekly,lag=k))/Var
}

plot(1:48,Gk,ylab='Gk',xlab='Lag',type='l')

#deseason
deseasoned <- weekly - weekly_components$seasonal
plot(deseasoned)

#deaseasoned stil has a trend. Need to difference
diff_ds <- diff(deseasoned,differences = 1)
adf.test(diff_ds,alternative = 'stationary')
plot(diff_ds,type='l',main='Time series differenced by 1',ylab='Differenced price')

acf(as.numeric(diff_ds))
pacf(as.numeric(diff_ds))

#fit some different arimas
fit <- arima(weekly,order=c(3,1,2),seasonal = list(order=c(1,0,0),period=52),method='CSS')
fit2 <- arima(weekly,order=c(2,1,2),seasonal = list(order=c(1,0,0),period=52),method='CSS')
fit3 <- arima(weekly,order=c(3,1,4),seasonal = list(order=c(1,1,0),period=52),method='CSS')
fit4 <- arima(weekly,order=c(0,1,2),seasonal = list(order=c(1,1,0),period=52),method='CSS')

#summary of fit
summary(fit);summary(fit2);summary(fit3);summary(fit4)

#diagnostic plots
tsdisplay(residuals(fit), main='(3,1,2)(1,0,0) Model Residuals',lag.max = 20)
tsdisplay(residuals(fit2), main='(2,1,2)(1,0,0) Model Residuals',lag.max=20)
tsdisplay(residuals(fit3), main='(3,1,4)(1,1,0) Model Residuals',lag.max=20)
tsdisplay(residuals(fit4), main='(0,1,2)(1,1,0) Model Residuals',lag.max=20)

#create a test set and training set
test <- window(as.numeric(weekly), start=183)
train <-window(as.numeric(weekly[1:182]))
train <- ts(train,frequency = 52)

#train models on subset of data
train_fit3 <- arima(train,order=c(0,1,2),seasonal = list(order=c(1,0,1),period=52),method='CSS')
train_fit4 <- arima(train,order=c(0,1,2),seasonal = list(order=c(1,1,0),period=52),method='CSS')
train_fit5 <- arima(train,order=c(0,1,1),seasonal = list(order=c(1,1,1),period=52),method='CSS')
train_fit6 <- arima(train,order=c(2,1,2),seasonal = list(order=c(1,0,1),period=52),method='CSS')
train_fit7 <- arima(train,order=c(2,1,1),seasonal = list(order=c(1,0,1),period=52),method='CSS')
holtWint <- HoltWinters(train,seasonal = 'additive')


#get forecats for holdout set
fcast_train3 <- forecast(train_fit3,h=20)
fcast_train4 <- forecast(train_fit4,h=20)
fcast_train5 <- forecast(train_fit5,h=20)
fcast_train6 <- forecast(train_fit6,h=20)
fcast_train7 <- forecast(train_fit7,h=20)
hwint_fcast <- forecast(holtWint,h=20)

par(mfrow=c(3,2))
#plot forecast and holdout data
plot(fcast_train3)
lines(weekly,col='red')


plot(fcast_train4)
lines(weekly,col='red')

plot(fcast_train5)
lines(weekly,col='red')

plot(fcast_train6)
lines(weekly,col='red')

plot(hwint_fcast)
lines(weekly,col='red')

plot(fcast_train7)
lines(weekly,col='red')

accuracy(fcast_train3,test);accuracy(fcast_train4,test);
accuracy(fcast_train5,test);accuracy(fcast_train6,test);
accuracy(fcast_train7,test);accuracy(hwint_fcast,test)


########### Feature Selection ###########

df$price = as.numeric(df$price)
drops = c("weekly_price", "square_feet", "monthly_price", "is_business_travel_ready")
df = df[, !(names(df) %in% drops)]

data.complete = df[complete.cases(df),]
data.complete$price = as.numeric(data.complete$price)
data.complete1 = data.matrix(data.complete)

str(data.complete)

factorvars = colnames(data.complete[,sapply(data.complete, is.factor)])
factorvars


#LASSO for feature selection since orginal dataframe had 96 explanatory variables

removevars = c(factorvars, "price", "latitude", "longitude")

y <- as.vector(data.complete$price)
x = data.complete[, !(names(data.complete) %in% removevars)]
# x = data.frame(matrix(unlist(x), nrow = nrow(x), byrow = T))
# x = as.data.frame(x)
xm = as.matrix(x)
# alpha = 1 - LASSO
lambdas <- exp( seq(-3, 10, length=50))
a <- glmnet(x=xm, y=y, lambda=rev(lambdas), 
            family='gaussian', alpha=1, intercept=TRUE)

plot(a, xvar='lambda', label=TRUE, lwd=6, cex.axis=1.5, cex.lab=1.2)
b = lars(x=xm, y=y, type="lasso", intercept=TRUE)
plot(b, lwd = 4)

##### ARIMAX model ########

#load data and clean it up
vars <-c('date',
         'price',
         'latitude',
         'longitude',
         'cleaning_fee',
         'review_scores_communication',
         'review_scores_accuracy',
         'calculated_host_listings_count',
         'availability_365',
         'review_scores_checkin',
         'reviews_per_month',
         'extra_people',
         'guests_included',
         'maximum_nights', 'last_review')


df2 <- df[,names(df) %in% vars]

#change to date time object
df2$last_review <- as.Date(df2$last_review, format="%Y-%m-%d")

#subset dataframe
df2 <- df2[df2$last_review > as.Date("2015-01-01"),]

#clean it up
df2 <- NaRV.omit(df2)
df2$cleaning_fee <- as.numeric(df2$cleaning_fee)
df2$extra_people <- as.numeric(df2$extra_people)
df2$price <- as.numeric(df2$price)
df3 <-  NaRV.omit(df2)
df3 <- as.xts(df3,order.by=as.Date(df3$last_review))
weeklydf <- apply.weekly(df3,mean)

n = length(weekly)
ntest = 20
ntrain = n-ntest

train.wk = weekly[1:ntrain]
test.wk  = weekly[(ntrain+1):n]

train.wk = ts(train.wk, frequency = 52)
test.wk = ts(test.wk, frequency = 52)

vtrain = weeklydf[1:ntrain,]
vtest  = weeklydf[(ntrain+1):n,]

fitarimax = arima(train.wk, order=c(2,1,2),
                  seasonal=list(order=c(1,0,1), period=52),
                  xreg=vtrain[,c(1:5, 7:8, 11:14)], method="CSS")

# acf(fitarimax$residuals)

p1 = predict(fitarimax, n.ahead=20, newxreg=vtest[,c(1:5, 7:8, 11:14)])

plot(p1$pred)

#train models on subset of data
train_fit3 <- arima(train.wk,order=c(0,1,2),seasonal = list(order=c(1,0,1),period=52), xreg=vtrain[,c(1:5, 7:8, 11:14)], method='CSS')
train_fit4 <- arima(train.wk,order=c(0,1,2),seasonal = list(order=c(1,1,0),period=52), xreg=vtrain[,c(1:5, 7:8, 11:14)], method='CSS')
train_fit5 <- arima(train.wk,order=c(0,1,1),seasonal = list(order=c(1,1,1),period=52), xreg=vtrain[,c(1:5, 7:8, 11:14)], method='CSS')
train_fit6 <- arima(train.wk,order=c(2,1,2),seasonal = list(order=c(1,0,1),period=52), xreg=vtrain[,c(1:5, 7:8, 11:14)], method='CSS')
train_fit7 <- arima(train.wk,order=c(2,1,1),seasonal = list(order=c(1,0,1),period=52), xreg=vtrain[,c(1:5, 7:8, 11:14)], method='CSS')


#get forecats for holdout set
fcast_train3 <- forecast(train_fit3, h=20, xreg=vtest[,c(1:5, 7:8, 11:14)])
fcast_train4 <- forecast(train_fit4, h=20, xreg=vtest[,c(1:5, 7:8, 11:14)])
fcast_train5 <- forecast(train_fit5, h=20, xreg=vtest[,c(1:5, 7:8, 11:14)])
fcast_train6 <- forecast(train_fit6, h=20, xreg=vtest[,c(1:5, 7:8, 11:14)])
fcast_train7 <- forecast(train_fit7, h=20, xreg=vtest[,c(1:5, 7:8, 11:14)])


par(mfrow=c(3,2))
#plot forecast and holdout data
plot(fcast_train3, main="Forecasts from ARIMAX(0,1,2)(1,0,1)[52]")
lines(weekly,col='red')


plot(fcast_train4, main="Forecasts from ARIMAX(0,1,2)(1,1,0)[52]")
lines(weekly,col='red')

plot(fcast_train5, main="Forecasts from ARIMAX(0,1,1)(1,1,1)[52]")
lines(weekly,col='red')

plot(fcast_train6, main="Forecasts from ARIMAX(2,1,2)(1,0,1)[52]")
lines(weekly,col='red')

plot(fcast_train7, main="Forecasts from ARIMAX(2,1,1)(1,0,1)[52]")
lines(weekly,col='red')

test.wk <- window(as.numeric(weekly), start=183)

accuracy(fcast_train3,test.wk);
accuracy(fcast_train4,test.wk);
accuracy(fcast_train5,test.wk);
accuracy(fcast_train6,test.wk);
accuracy(fcast_train7,test.wk)

##### DONE ########
