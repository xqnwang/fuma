`fuma`
========

Forecast uncertainty based on model averaging.

The R package `fuma` provides implementations of the uncertainty estimation of feature-based time series forecasts, see our [paper](https://arxiv.org/abs/1908.02891) for the details.

Installation
------------

You can install the package `fuma` from [GitHub Repository](https://github.com/xqnwang/fuma) with:

``` r
devtools::install_github("xqnwang/fuma")
```

Usage
-----

This part explains how to reproduce the results for our working paper. The feature-based framework for the uncertainty estimation is applied in this part to obtain the forecasts and prediction intervals of M3 dataset for the confidence level $95\%$.

### The reference dataset and test dataset

``` r
library(fuma)
library(Mcomp)
library(magrittr)

# summarize the distribution of sample size on the M3 dataset
yearly_M3 <- Filter(function(l) l$period == "YEARLY", M3)
quarterly_M3 <- Filter(function(l) l$period == "QUARTERLY", M3)
monthly_M3 <- Filter(function(l) l$period == "MONTHLY", M3)
l_y <- sapply(yearly_M3, function(lentry)lentry$n)
l_q <- sapply(quarterly_M3, function(lentry)lentry$n)
l_m <- sapply(monthly_M3, function(lentry)lentry$n)

# generate reference dataset by GRATIS
set.seed(2019-03-11-1)
x_y <- ts_generate(n.ts = 10000, freq = 1, length = l_y, h = 6)
set.seed(2019-03-11-2)
x_q <- ts_generate(n.ts = 10000, freq = 4, length = l_q, h = 8)
set.seed(2019-03-11-3)
x_m <- ts_generate(n.ts = 10000, freq = 12, length = l_m, h = 18)
train <- append(x_y, x_q) %>% append(x_m)

# test dataset (M3 dataset)
test <- Filter(function(l) l$period != "OTHER", M3)
```

### Train benchmark models on reference dataset

``` r
# train benchmark models and obtain forecasts and prediction intervals
train_forec <- ts_forec(dataset = train, methods_list(), level = c(80, 85, 90, 95, 99), parallel = TRUE)

# evaluate forecasting performance
train_scores <- calc_scores(train_forec, parallel = TRUE) 

# extract the MSIS scores
level <- 95
train_msis <- extract_msis(train_scores, level = level) 
```

### Extract time series features from reference dataset

``` r
# scale the time series
train_scaledx <- lapply(train, function (lentry){
  mu <- mean(lentry$x)
  sigma <- sd(lentry$x)
  x <- (lentry$x - mu)/sigma
  lentry[names(lentry) == "x"] <- list(x)
  lentry
})

# calculate features
train_feat0 <- M4metalearning::THA_features(train_scaledx)
train_feat <- t(sapply(train_feat0, function (lentry) {
  seriesdata <- c(as.numeric(lentry$features))
  names(seriesdata) <- c(names(lentry$features))
  seriesdata
}))
```

### Linking features with interval forecasting performance

The relationship between features and interval forecasting accuracy is captured by generalized additive model in this paper because of its properties of interpretability, flexibility, automation and regularization.

``` r
# combine features and MSIS values
train_fm <- `names<-` (append(list(train_feat), list(train_msis)), 
                       c("feat", "msis"))

X <- as.data.frame(train_fm$feat)
Y <- as.data.frame(train_fm$msis)
X$seasonal_period_q <- ifelse(X$seasonal_period == 4, 1, 0)
X$seasonal_period_m <- ifelse(X$seasonal_period == 12, 1, 0)
X <- subset(X, select = -seasonal_period)

# GAM model training
gam_fit <- gam.fun(X, Y, LogY = TRUE, k = 10, parallel = TRUE)
```

### Optimal threshold search

``` r
# get fitted value
gam_fitvalue <- matrix(NA, nrow = nrow(X), ncol = ncol(Y)) %>% data.frame()
rownames(gam_fitvalue) <- rownames(X)
colnames(gam_fitvalue) <- colnames(Y)
for (i in 1:length(gam_fit)){
  gam_fitvalue[, i] <- gam_fit[[i]]$fitted.values
}

# adjusted softmax transformation
gam_prob <- softmax.fun(gam_fitvalue)

# search optimal threshold
init_ratio <- seq(0.1, 1, 0.1)
gam.mean_threshold <- choose.threshold(dataset = train_forec, fitProb = gam_prob, 
                                       ratio = init_ratio, combine="mean", 
                                       level = level, parallel = TRUE) 
gam.weighted_threshold <- choose.threshold(dataset = train_forec, fitProb = gam_prob, 
                                           ratio = init_ratio, combine="weighted", 
                                           level = level, parallel = TRUE) 

# optimal threshold
gam.mean_ratio <- c(gam.mean_threshold$threshold_y, 
                    gam.mean_threshold$threshold_q, 
                    gam.mean_threshold$threshold_m)
gam.weighted_ratio <- c(gam.weighted_threshold$threshold_y, 
                        gam.weighted_threshold$threshold_q, 
                        gam.weighted_threshold$threshold_m)
```

### Forecast test dataset

``` r
# extract features from test dataset
test_scaledx <- lapply(test, function (lentry){
  mu <- mean(lentry$x)
  sigma <- sd(lentry$x)
  x <- (lentry$x - mu)/sigma
  lentry[names(lentry) == "x"] <- list(x)
  lentry
})
test_feat0 <- THA_features(test_scaledx)
test_feat <- t(sapply(test_feat0, function (lentry) {
  seriesdata <- c(as.numeric(lentry$features))
  names(seriesdata) <- c(names(lentry$features))
  seriesdata
}))

# dummy out categorical features
Xtest <- as.data.frame(test_feat)
Xtest$seasonal_period_q <- ifelse(Xtest$seasonal_period == 4, 1, 0)
Xtest$seasonal_period_m <- ifelse(Xtest$seasonal_period == 12, 1, 0)
Xtest <- subset(Xtest, select = -seasonal_period)

# predicted values of GAMs and weight assignment
gam_pre <- matrix(NA, nrow = nrow(Xtest), ncol = ncol(Y)) %>% data.frame()
rownames(gam_pre) <- rownames(Xtest)
colnames(gam_pre) <- colnames(Y)
for (i in 1:length(gam_fit)){
  gam_pre[, i] <- predict(gam_fit[[i]], newdata = Xtest)
}
gam_preprob <- softmax.fun(gam_pre)

# forecasts of benchmark methods
test_forec <- ts_forec(dataset = test, methods_list(), level = 95, parallel = TRUE)

# OurMethod(mean)
fuma_forec <- comb_forec(test_forec, weightprob = gam_preprob,
                         Threshold = gam.mean_ratio, level = level,
                         combine = "mean", methodname = "OurMethod(mean)",
                         show.methods = FALSE, parallel = TRUE)

# OurMethod(weighted)
fuma_forec <- comb_forec(fuma_forec, weightprob = gam_preprob,
                         Threshold = gam.weighted_ratio, level = level,
                         combine = "weighted", methodname = "OurMethod(weighted)",
                         show.methods = FALSE, parallel = TRUE)

# OurMethod(allweighted)
all_ratio <- c(-1, -1, -1) 
fuma_forec <- comb_forec(fuma_forec, weightprob = gam_preprob,
                         Threshold = all_ratio, level = level,
                         combine = "weighted", methodname = "OurMethod(allweighted)",
                         show.methods = FALSE, parallel = TRUE)
```

References
----------

- Wang, X., Kang, Y., Petropoulos, F., & Li, F. (2019). Que ser\'a ser\'a? The uncertainty estimation of feature-based time series forecasts. [Working paper on arXiv](https://arxiv.org/abs/1908.02891).


