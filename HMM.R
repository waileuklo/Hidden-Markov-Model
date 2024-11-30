library(depmixS4)
library(nnet)
library(MASS)
library(Rsolnp)
library(nlme)
library(TTR)
library(ggplot2)
library(reshape2)
library(xts)
library(quantmod)
library(Quandl)
library(fHMM)
library(tidyverse)
library(fUnitRoots)

################################### Tutorial #######################################
getSymbols("^TWII", src="yahoo", from="1900-01-01", to="2020-01-03")
chartSeries(TWII, theme="black")
TWII_subset = window(TWII, start=as.Date("2018-01-01"), end=as.Date("2020-01-13"))
TWII_train = cbind(TWII_subset$TWII.Close-TWII_subset$TWII.Open)
TWII_train
# 
# mod = depmix(TWII.Close~1, family=gaussian(),
#              nstates=5, data=TWII_train)
# mod@response
# mod@transition
# mod@prior
# mod@ntimes
# set.seed(1)
# fm2 = fit(mod, verbose=FALSE)
# fm2@transition
# 
# probs = posterior(fm2)
# head(probs)
# probs$state
# TWII_Predict = cbind(TWII_subset$TWII.Close, probs$state)
# 
# chartSeries(TWII_Predict[,1])
# addTA(TWII_Predict[TWII_Predict[,2]==1,1],on=1,type="p",col=5,pch=25)
# addTA(TWII_Predict[TWII_Predict[,2]==2,1],on=1,type="p",col=6,pch=24)
# addTA(TWII_Predict[TWII_Predict[,2]==3,1],on=1,type="p",col=7,pch=23)
# addTA(TWII_Predict[TWII_Predict[,2]==4,1],on=1,type="p",col=8,pch=22)
# addTA(TWII_Predict[TWII_Predict[,2]==5,1],on=1,type="p",col=10,pch=21)
# 
# head(TWII_Predict)


################################### GARCH #######################################
setwd('C:/Users/johnl/My Drive/Documents/HKUST/BDT/Courses/MSBD 5006/Project/data')

data = read.csv('gold_dec24(GC=F)_1d.csv', header=TRUE)
data$Date = as.Date(data$Date, format="%Y-%m-%d %H:%M:%S")

# data = xts(data$Close, data$Date)
# colnames(data) = c("Close")
# head(data)
# data_train = window(data, start=as.Date("2004-11-01"), end=as.Date("2024-04-30"))
# chartSeries(data_train, theme="black")

data_train = data.frame(data$Date, data$Close)
colnames(data_train) = c("date", "price")
length(data_train$price)
unitrootTest(data_train$price,lags=1,type=c("c"))

par(mfrow=c(1,1))
plot(data_train$price,xlab=" ",ylab=" ",main="Daily Gold Close Price",type="l",col="black", lwd=1)

# simple return
simple_return = data_train$price[2:length(data_train$price)] / data_train$price
unitrootTest(simple_return,lags=1,type=c("c"))
# The p-value is smaller than 0.05, so reject
# null and conclude that the simple return
# is stationary


# log return
log_return = diff(log(data_train$price))
# percentage_log_return = percentage_log_return[2:length(percentage_log_return)]
plot(log_return,xlab=" ",ylab=" ",main="Daily log returns for Gold",type="l",col="black", lwd=1)

# percentage_log_return_all = diff(log(data))*100
# percentage_log_return_all = percentage_log_return_all[2:length(percentage_log_return_all)]
# plot(percentage_log_return_all,xlab=" ",ylab=" ",main="Daily percentage log returns for Gold",type="l",col="red")

# There is serial correlation in the log return
acf(log_return,20,main="20 lags",col="red")
pacf(log_return,20,main="20 lags",col="red")
Box.test(log_return,10,type="Ljung")

# There is ARCH effect in the log return
at=log_return - mean(log_return)
Box.test(at^2,10,type="Ljung")

library(rugarch)
library(fGarch)
spec1=ugarchspec(variance.model=list(model="iGARCH"),
                 mean.model=list(armaOrder=c(0,0),include.mean=TRUE),
                 distribution.model="norm")
mm=ugarchfit(spec=spec1,data=log_return)
# mm=garchFit(~arma(2,1)+garch(1,1),data=log_return,include.mean=FALSE,trace=F)
mm  ### see output

res=residuals(mm,standardize=T)
Box.test(res,10,type="Ljung")
Box.test(res,20,type="Ljung")
Box.test(res^2,10,type="Ljung")

pred = ugarchforecast(mm, n.ahead=4)
pred_log_return = fitted(pred)
pred_se = sigma(pred)
pred_ub = (fitted(pred) + 1.96 * pred_se)[1]
pred_lb = (fitted(pred) - 1.96 * pred_se)[1]


fore = predict(mm,n.ahead=7)
fore

L1=append(coredata(log_return)[4897],fore[1,1]-1.96*fore[1,3])
L1
U1=append(coredata(log_return)[4897],fore[1,1]+1.96*fore[1,3])
U1

par(mfrow=c(1,1))
plot(1:21,log_return[4884:4904],type="o",ylab="",xlab="",ylim=c(-0.05,0.05),main="Forecasting by iGARCH")
lines(14:21,append(coredata(log_return)[4897],fore[,1]),type="o",col="red")
lines(14:15, U1,type="l",col="blue")
lines(14:15, L1,type="l",col="blue")
#points(temp.date[2:7],p,type="o")
legend(x="topleft",c("True returns","prediction"),lty=c(1,1),pch=c(1,1),col=c("black","red"))

# Rolling prediction
bundle = function(x){
  # the log return from day (x-1) determines the price on day x
  log_return_subset = log_return[(x-1):(x+89)]
  # print(log_return_subset)
  
  # model = garchFit(~arma(2,1)+garch(1,1),data=log_return_subset,include.mean=FALSE,trace=F)
  # pred = predict(model,n.ahead=1)
  # pred_log_return = pred[1, 1]
  # pred_se = pred[1, 3]
  # pred_ub = pred_log_return + 1.96 * pred_se
  # pred_lb = pred_log_return - 1.96 * pred_se
  
  spec = ugarchspec(variance.model=list(model="iGARCH"),
                    mean.model=list(armaOrder=c(0,0),include.mean = TRUE),
                    distribution.model="norm")
  model = ugarchfit(spec=spec,data=log_return_subset)
  pred = ugarchforecast(model, n.ahead=1)
  pred_log_return = fitted(pred)
  pred_se = sigma(pred)
  pred_ub = (fitted(pred) + 1.96 * pred_se)[1]
  pred_lb = (fitted(pred) - 1.96 * pred_se)[1]
  
  pred_E = exp(pred_log_return + pred_se^2/2)
  pred_U = exp(pred_ub)
  pred_L = exp(pred_lb)

  last_price = data_train$price[x + 90]
  next_price = data_train$price[x + 91]
  pred_E = last_price * pred_E
  pred_U = last_price * pred_U
  pred_L = last_price * pred_L
  c(next_price, pred_E, pred_U, pred_L)
  
  ## compare true log return with predicted log return
  # next_log_return = log_return[x + 90]
  # c(next_log_return, pred_log_return, pred_ub, pred_lb)
  
  ## compare true gross return with predicted gross return
  # next_gross_return = data_train$price[x + 91] / data_train$price[x + 90]
  # c(next_gross_return, pred_E, pred_U, pred_L)
  
}

start = 4816
end = 4934
result = t(sapply(start:end, bundle))
colnames(result) = c("price", "prediction", "upper bound", "lower bound")
result = data.frame(result)
mape = mean((abs(result[,1] - result[,2])) / result[,1] * 100)
mape
period = as.Date(data_gold$dates[(start+91):(end+91)])

# period_to_plot = 1:length(period)
period_to_plot = 95:100
forecast = data.frame(
  Time=period[period_to_plot],
  Real=result$price[period_to_plot],
  Pred=result$prediction[period_to_plot],
  Upper=result$upper.bound[period_to_plot],
  Lower=result$lower.bound[period_to_plot]
)
ggplot(forecast, aes(x=Time, y=Pred)) +
  ggtitle("Prediction of Gold Price in Year 2024") +
  geom_ribbon(aes(x=Time, ymin=Upper, ymax=Lower), fill = "grey", alpha = 0.9) +
  geom_line(aes(y=Real, colour="historic price")) +
  geom_line(aes(y=Pred, colour="prediction"), alpha=0.5) +
  scale_color_manual(name="legend", values=c("historic price"="black", "prediction"="blue"))




################################### depmixS4 ###################################
# Hidden Markov Model for gold price
hmm = depmix(Close~1, family=gaussian(),
             nstates=4, data=data_train)
hmm
set.seed(1)
fhmm = fit(hmm, verbose=FALSE)
fhmm@transition

probs = posterior(fhmm)
head(probs)
probs$state
data_Predict = cbind(data_train$Close, probs$state)

chartSeries(data_Predict[,1])
addTA(data_Predict[data_Predict[,2]==1,1],on=1,type="p",col=5,pch=25)
addTA(data_Predict[data_Predict[,2]==2,1],on=1,type="p",col=6,pch=24)
addTA(data_Predict[data_Predict[,2]==3,1],on=1,type="p",col=7,pch=23)
addTA(data_Predict[data_Predict[,2]==4,1],on=1,type="p",col=8,pch=22)

# extract the state-transition matrix
transition_mat <- rbind(getpars(getmodel(fhmm,"transition",1)),getpars(getmodel(fhmm,"transition",2)))
transition_mat

# extract the probability of the states at the final time point in the data (t=T)
# this will act as a "prior" to compute the forecasted state distributions
prior_vec <- as.numeric(posterior(fhmm)[1000,-1])
prior_vec

# state-wise predictions for the observed variables
pred_r_by_state <- c(getpars(getmodel(fhmm,"response",1))[1],
                     getpars(getmodel(fhmm,"response",2))[1])
pred_r_by_state

# for T + 1
# the forecasted state distribution is (prior_vec %*% transition_mat)
# so hence the prediction of the observed variable is
sum(pred_r_by_state * (transition_mat %*% prior_vec))

# for T + 2
# the forecasted state distribution is (prior_vec %*% transition_mat %*% transition_mat)
# so hence the prediction of the observed variable is
sum(pred_r_by_state * (prior_vec %*% transition_mat %*% transition_mat))

# for T + 3
sum(pred_r_by_state * (prior_vec %*% transition_mat %*% transition_mat %*% transition_mat))


################################### fHMM #######################################
dev.off()
setwd('C:/Users/johnl/My Drive/Documents/HKUST/BDT/Courses/MSBD 5006/Project/data')

# Load data
controls = list(states=3, sdds="t",
                data = list(file='gold_dec24(GC=F)_1d.csv', date_column="Date", data_column="Close",
                            logreturns=TRUE, from="2004-11-01", to="2024-10-30"),
                fit=list(runs=100))
controls <- set_controls(controls)
class(controls)
data_gold <- prepare_data(controls)
summary(data_gold)
data_gold_ts = data.frame(data_gold$dates, data_gold$time_series)
colnames(data_gold_ts) = c("date", "price")

plot(data_gold)

# Fit model
model_3n <- fit_model(data_gold, seed=42)
coef(model_3n, alpha=0.05)
par(mfrow=c(1,2))
plot(model_3n, plot_type=c("ll", "sdds"), colors=c('firebrick2', 'gold1', 'cadetblue1'))
summary(model_3n)

# Visualize hidden states
model_3n = decode_states(model_3n)
table(model_3n$decoding)
events <- fHMM_events(list(
  dates = c("2008-09-15"),
  labels = c("Bankruptcy of Lehman Brothers")))
par(mfrow=c(1,1))
plot(model_3n, colors=c('firebrick2', 'gold1', 'cadetblue1'),
     events=events,
     from="2007-11-01", to="2010-01-30")

# Model checking
model_3n = compute_residuals(model_3n)
plot(model_3n, plot_type = "pr")
res = residuals(model_3n)
res
Box.test(res,10,type="Ljung")
tseries::jarque.bera.test(res)


# Rolling prediction
start = "2004-11-01"
bundle = function(x){
  from = data_gold$dates[x]
  to = data_gold$dates[x+90]
  controls = list(states=3, sdds="t",
                  data = list(file='gold_dec24(GC=F)_1d.csv', date_column="Date", data_column="Close",
                              logreturns=TRUE, from=from, to=to),
                  fit=list(runs=100))
  controls = set_controls(controls)
  data_subset = prepare_data(controls)
  model = fit_model(data_subset, seed=42, verbose=FALSE)
  model = decode_states(model)
  
  pred = predict(model, ahead = 1)
  pred_log_return = pred$data[,2]
  pred_lb = pred$data[,1]
  pred_ub = pred$data[,3]
  pred_se = (pred_ub - pred_log_return) / 1.96
  
  pred_E = exp(pred_log_return + pred_se^2/2)
  pred_U = exp(pred_ub)
  pred_L = exp(pred_lb)
  
  len = length(data_subset$time_series)
  last_price = data_subset$time_series[len]
  # next_date = format(as.Date(start) + x + 91, "%Y-%m-%d")
  # next_price = data_gold_ts[data_gold_ts$date == next_date,]$price
  next_price = data_gold_ts[x+91,]$price
  pred_E = last_price * pred_E
  pred_U = last_price * pred_U
  pred_L = last_price * pred_L
  c(next_price, pred_E, pred_U, pred_L)
}

start = 4816
end = 4934
sequence = seq(start, end, by=1)
result = t(sapply(sequence, bundle))
colnames(result) = c("price", "prediction", "upper bound", "lower bound")
result = data.frame(result)
mape = mean((abs(result[,1] - result[,2])) / result[,1] * 100)
mape
period = as.Date(data_gold$dates[(start+91):(end+91)])

period_to_plot = 1:length(period)
# period_to_plot = 95:100
forecast = data.frame(
  Time=period[period_to_plot],
  Real=result$price[period_to_plot],
  Pred=result$prediction[period_to_plot],
  Upper=result$upper.bound[period_to_plot],
  Lower=result$lower.bound[period_to_plot]
)
ggplot(forecast, aes(x=Time, y=Pred)) +
  ggtitle("Prediction of Gold Price in Year 2024") +
  geom_ribbon(aes(x=Time, ymin=Upper, ymax=Lower), fill = "grey", alpha = 0.9) +
  geom_line(aes(y=Real, colour="historic price")) +
  geom_line(aes(y=Pred, colour="prediction"), alpha=0.5) +
  scale_color_manual(name="legend", values=c("historic price"="black", "prediction"="blue"))


# Prediction test
pred = predict(model_3n, ahead = 4)
pred_log_return = pred$data[,2]
pred_lb = pred$data[,1]
pred_ub = pred$data[,3]
pred_se = (pred_ub - pred_log_return) / 1.96

pred_E = exp(pred_log_return + pred_se^2/2)
pred_U = exp(pred_ub)
pred_L = exp(pred_lb)



controls = list(states=3, sdds="t",
                data = list(file='gold_dec24(GC=F)_1d.csv', date_column="Date", data_column="Close",
                            logreturns=FALSE, from="2019-03-01", to="2020-11-06"),
                fit=list(runs=100))
controls <- set_controls(controls)
data_gold_extend <- prepare_data(controls)
data_gold_extend = data_gold_extend$data
return = diff(data_gold_extend)
len = length(data_gold$data)
forecast = append(return[len], pred_E)
U = append(return[len], pred_U)
L = append(return[len], pred_L)

par(mfrow=c(1,1))
plot(1:16,return[(len-11):(len+4)],type="o",ylab="",xlab="",main="4-step Forecast")
lines(12:16,forecast,type="o",col="red")
lines(12:16, U,type="l",col="blue")
lines(12:16, L,type="l",col="blue")
legend(x="bottomleft",c("True return","prediction"),lty=c(1,1),pch=c(1,1),col=c("black","red"))

