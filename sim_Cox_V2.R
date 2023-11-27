# play weibull distribution
N=100
shape <- 8
scale <- 5
t <- scale * (-log(runif(N)))^(1/shape)
hist(t, probability = T, xlim=c(0,20))
lines(density(t))
u <- seq(0, 20, 0.1)
lines(u, dweibull(u, shape=8, scale=5), col="red")
lines(u, dexp(u, 0.06), col="blue")

library(MASS)
library(glmnet)
library(xrnet)
library(ggplot2)
# Rcpp::sourceCpp('cindex.cpp')
# Rcpp::sourceCpp('cdl_cox_rcpp.cpp')

# Simulation cox model with external information
N <- 100
p <- 200
q <- 10
shape <- 8
scale <- 5
nonzero_a <- 0.2
effect_size <- 0.2
set.seed(201909)
# generate z, alpha
z <- rbinom(p*q, size = 1, prob = 0.1)
z <- matrix(z, nrow = p, ncol = q)
mode(z) <- "double"
a <- rep(0, q)
a[c(1:floor(nonzero_a*q))] <- effect_size
# generate beta
snr_z <- 2
var_za <- var(drop(z %*% a))  # empirical variance
var_b <- var_za / snr_z
b <- drop(z %*% a) + sqrt(var_b)*rnorm(p)
# generate x
corr_x <- 0.5
sigma <- matrix(NA, nrow = p, ncol = p)
for (j in 1:p) {
    for (i in 1:p) {
        sigma[i, j] <- corr_x^abs(i - j)
    }
}
x <- mvrnorm((1e4+N), rep(0, p), sigma)

start <- Sys.time()
# generate survival time, fix true c index at 0.8
for (stdev in seq(0.01, 10, 0.01)) {
    t <- scale * (-log(runif(1e4+N)) * exp(-drop(x%*%b)))^(1/shape) + rnorm(n=(1e4+N), sd=stdev)
    t[t<=0] <- min(t[t>0])
    t[t>20] <- 20    # follow up time is 20
    c <- rexp((1e4+N), 0.06)     # produce ratio of event to censor at around 3/1
    c[c>20] <- 20
    time <- pmin(t,c)
    status <- ifelse(t<c, 1, 0)
    table(status)
    y <- cbind(time, status)
    Cindex(drop(x%*%b), y)
    if ( abs(Cindex(drop(x%*%b), y)-0.8) < 0.001 ) break
}
end <- Sys.time()
end-start
print(stdev)
print(Cindex(drop(x%*%b), y))

x_val <- x[1:N, ]
y_val <- y[1:N, ]
x_test <- x[-(1:N), ]
y_test <- y[-(1:N), ]


x <- mvrnorm(N, rep(0, p), sigma)
t <- scale * (-log(runif(N)) * exp(-drop(x%*%b)))^(1/shape) + rnorm(N, sd = stdev)
t[t<=0] <- min(t[t>0])
t[t>20] <- 20    # follow up time is 20
c <- rexp(N, 0.06)     # produce ratio of event to censor at around 3/1
c[c>20] <- 20
time <- pmin(t,c)
status <- ifelse(t<c, 1, 0)
table(status)
y <- cbind(time, status)

ordered <- order(y[,1])
y <- y[ordered, ]
x <- as.matrix(x[ordered, ])

xm <- colMeans(x)
x_norm <- sweep(x, 2, xm, "-")
xs <- drop(sqrt(crossprod(rep(1/N,N), x_norm^2)))
x_norm <- sweep(x_norm, 2, xs, "/")

fit_cdl = cdlcoxRcpp(y, x_norm, z, standardize = F)
fit_xrnet = xrnet(x_norm, y, z, family = "cox", 
                  penalty_main = define_ridge(), penalty_external = define_lasso(penalty_ratio = 0.01),
                  standardize = c(F,F))

rbind(fit_cdl$lambda1, fit_xrnet$penalty)
rbind(fit_cdl$lambda2, fit_xrnet$penalty_ext)

rbind(fit_cdl$beta[1:10,3], fit_xrnet$betas[1:10,3,1])
rbind(fit_cdl$beta[201:210,3], fit_xrnet$alphas[,3,1])

rbind(fit_cdl$beta[1:10,43], fit_xrnet$betas[1:10,3,3])
rbind(fit_cdl$beta[201:210,43], fit_xrnet$alphas[,3,3])

rbind(fit_cdl$beta[1:10,190], fit_xrnet$betas[1:10,10,10])
rbind(fit_cdl$beta[201:210,190], fit_xrnet$alphas[,10,10])

rbind(fit_cdl$beta[1:10,400], fit_xrnet$betas[1:10,20,20])
rbind(fit_cdl$beta[201:210,400], fit_xrnet$alphas[,20,20])
#########################################################################################################################
nSim <- 100
cidx_test <- numeric()
coefs <- mat.or.vec(nSim, p)
cidx_test_noext <- numeric()
coefs_noext <- mat.or.vec(nSim, p)
tpr <- numeric()
fpr <- numeric()
alphas <- mat.or.vec(nSim, q)
for (j in 1:nSim) {
    # xrnet
    x <- mvrnorm(N, rep(0, p), sigma)
    t <- scale * (-log(runif(N)) * exp(-drop(x%*%b)))^(1/shape) + rnorm(N, sd = stdev)
    t[t<=0] <- min(t[t>0])
    t[t>20] <- 20    # follow up time is 20
    c <- rexp(N, 0.06)     # produce ratio of event to censor at around 3/1
    c[c>20] <- 20
    time <- pmin(t,c)
    status <- ifelse(t<c, 1, 0)
    table(status)
    y <- cbind(time, status)
    
    ordered <- order(y[,1])
    y <- y[ordered, ]
    x <- as.matrix(x[ordered, ])
    
    fit <- xrnet(x, y, z, family="cox", 
                 penalty_main=define_ridge(penalty_ratio=0.01), 
                 penalty_external=define_lasso(penalty_ratio=0.01))
    beta <- fit$betas
    alpha <- fit$alphas
    cidx_val <- apply(beta, c(2,3), function(l) Cindex(drop(x_val%*%l), y_val))
    wh <- drop(which(cidx_val==max(cidx_val), arr.ind = T))
    coefs[j, ] <- beta[, wh[1], wh[2]]
    alphas[j, ] <- alpha[, wh[1], wh[2]]
    tpr[j] <- sum(alphas[j, 1:floor(nonzero_a*q)]!=0) / floor(nonzero_a*q)
    fpr[j] <- sum(alphas[j, -(1:floor(nonzero_a*q))]!=0) / (q - floor(nonzero_a*q))
    cidx_test[j] <- Cindex(drop(x_test%*%beta[, wh[1], wh[2]]), y_test)

    # glmnet
    fit_glmnet <- glmnet(x, y, alpha = 0, family = "cox")
    beta_noext <- as.matrix(fit_glmnet$beta)
    cidx_val_noext <- apply(beta_noext, 2, function(l) Cindex(drop(x_val%*%l), y_val))
    wh_noext <- which.max(cidx_val_noext)
    coefs_noext[j, ] <- beta_noext[, wh_noext]
    cidx_test_noext[j] <- Cindex(drop(x_test%*%beta_noext[, wh_noext]), y_test)
}

boxplot(cidx_test, cidx_test_noext,
        col = c("blue", "cyan"), ylim = c(0.6, 0.8),
        names = c("RidgeLasso_high", "ridge_high"),
        main = "high external signal")
abline(h=0.8, lty=2)
#########################################################################################################################

#########################################################################################################################
source("fast_regularized_cox.R")
rm(cdl_Cox)

ThemeMain <- theme(axis.title.x = element_text(size=18),
                   axis.title.y = element_text(size=18),
                   axis.text.x = element_text(size=16),
                   axis.text.y = element_text(size=16),
                   legend.text = element_text(size=18),
                   legend.title = element_text(size=17),
                   plot.title = element_text(size=20, face='bold', hjust=0.5),
                   )

# base case: N = 100, p = 200, q = 50, true_cindex = 0.8, snr_z = 2, cov_x = "auto", alpha1 = 0, alpha2 = 1
# 1: high external signal vs. low, ridge-lasso
set.seed(1)
start <- Sys.time()
zeroext <- simCoxExt(N = 100, p = 200, q = 50, true_cindex = 0.8, snr_z = 5e-2, nSim = 100)
base <- simCoxExt(N = 100, p = 200, q = 50, true_cindex = 0.8, snr_z = 2, nSim = 100)
medext <- simCoxExt(N = 100, p = 200, q = 50, true_cindex = 0.8, snr_z = 1, nSim = 100)
lowext <- simCoxExt(N = 100, p = 200, q = 50, true_cindex = 0.8, snr_z = 0.1, nSim = 100)
end <- Sys.time()
end-start
a <- data.frame(cindex = base$cidx_test, signal=2, Method="Ridge_Lasso")
b <- data.frame(cindex = base$cidx_test_noext, signal=2, Method="Ridge")
c <- data.frame(cindex = medext$cidx_test, signal = 1, Method="Ridge_Lasso")
d <- data.frame(cindex = medext$cidx_test_noext, signal=1, Method="Ridge")
e <- data.frame(cindex = lowext$cidx_test, signal=0.1, Method="Ridge_Lasso")
f <- data.frame(cindex = lowext$cidx_test_noext, signal=0.1, Method="Ridge")
g <- data.frame(cindex = zeroext$cidx_test, signal=0, Method="Ridge_Lasso")
h <- data.frame(cindex = zeroext$cidx_test_noext, signal=0, Method="Ridge")
dat_plot <- rbind(a,b,c,d,e,f,g,h)
dat_plot$signal <- as.factor(dat_plot$signal)
cidx1 <- ggplot(data = dat_plot, mapping = aes(x = signal, y = cindex, fill = Method), position = "dodge") +
         geom_boxplot(width = 0.5) + 
         geom_hline(aes(yintercept = 0.8), linetype = "dashed") +
         scale_fill_manual(values = c("gray30", "gray80")) +
         xlab("meta-feature signal-noise ratio, SNR") +
         ylab("test C-index")
cidx1 <- cidx1 + ggtitle("Experiment 1 \n Prediction Performance") + ThemeMain; cidx1

a <- data.frame(tpr = base$tpr, signal = 2)
b <- data.frame(tpr = medext$tpr, signal = 1)
c <- data.frame(tpr = lowext$tpr, signal = 0.1)
d <- data.frame(tpr = zeroext$tpr, signal = 0)
dat_plot <- rbind(a,b,c,d)
dat_plot$signal <- as.factor(dat_plot$signal)
tpr1 <- ggplot(data = dat_plot, mapping = aes(x=signal, y=tpr, fill=signal), position = "dodge") +
        geom_boxplot(width = 0.5) + 
        scale_fill_manual(values = c("gray100", "gray60", "gray30", "gray10")) +
        xlab("meta-feature signal-noise ratio, SNR") +
        ylab("true postive rate") + 
        ggtitle("Experiment 1 \n Meta-Feature Selection True Positive Rate") + ThemeMain; tpr1

a <- data.frame(fpr = base$fpr, signal = 2)
b <- data.frame(fpr = medext$fpr, signal = 1)
c <- data.frame(fpr = lowext$fpr, signal = 0.1)
d <- data.frame(fpr = zeroext$fpr, signal = 0)
dat_plot <- rbind(a,b,c,d)
dat_plot$signal <- as.factor(dat_plot$signal)
fpr1 <- ggplot(data = dat_plot, mapping = aes(x=signal, y=fpr, fill=signal), position = "dodge") +
        geom_boxplot(width = 0.5) + 
        scale_fill_manual(values = c("gray100", "gray60", "gray30", "gray10")) +
        xlab("meta-feature signal-noise ratio, SNR") +
        ylab("false postive rate") + 
        ggtitle("Experiment 1 \n Meta-Feature Selection False Positive Rate") + ThemeMain; fpr1

a <- data.frame(acc = base$acc, signal = 2)
b <- data.frame(acc = medext$acc, signal = 1)
c <- data.frame(acc = lowext$acc, signal = 0.1)
d <- data.frame(acc = zeroext$acc, signal = 0)
dat_plot <- rbind(a,b,c,d)
dat_plot$signal <- as.factor(dat_plot$signal)
acc1 <- ggplot(data = dat_plot, mapping = aes(x=signal, y=acc, fill=signal), position = "dodge") +
    geom_boxplot(width = 0.5) + 
    scale_fill_manual(values = c("gray100", "gray60", "gray30", "gray10")) +
    xlab("meta-feature signal-noise ratio, SNR") +
    ylab("accuracy") + 
    ggtitle("Experiment 1 \n Meta-Feature Selection Accuracy") + ThemeMain; acc1
#########################################################################################################################
# 2: large sample size vs. small, ridge-lasso, high signal
set.seed(1)
mSample <- simCoxExt(N = 200, p = 200, q = 50, true_cindex = 0.8, snr_z = 2, nSim = 100)
lSample <- simCoxExt(N = 500, p = 200, q = 50, true_cindex = 0.8, snr_z = 2, nSim = 100)
a <- data.frame(cindex = base$cidx_test, sample.size="100", Method="Ridge_Lasso")
b <- data.frame(cindex = base$cidx_test_noext, sample.size="100", Method="Ridge")
c <- data.frame(cindex = mSample$cidx_test, sample.size="200", Method="Ridge_Lasso")
d <- data.frame(cindex = mSample$cidx_test_noext, sample.size="200", Method="Ridge")
e <- data.frame(cindex = lSample$cidx_test, sample.size="500", Method="Ridge_Lasso")
f <- data.frame(cindex = lSample$cidx_test_noext, sample.size="500", Method="Ridge")
dat_plot <- rbind(a,b,c,d,e,f)
cidx2 <- ggplot(data = dat_plot, mapping = aes(x = sample.size, y = cindex, fill = Method), position = "dodge") +
    geom_boxplot(width = 0.6) + 
    geom_hline(aes(yintercept = 0.8), linetype = "dashed") +
    scale_fill_manual(values = c("gray30", "gray80")) +    
    xlab("sample size, N") +
    ylab("test C-index") 
cidx2 <- cidx2 + ggtitle("Experiment 2 \n Prediction Performance") + ThemeMain; cidx2

a <- data.frame(tpr = base$tpr, sample.size = "100")
b <- data.frame(tpr = mSample$tpr, sample.size = "200")
c <- data.frame(tpr = lSample$tpr, sample.size = "500")
dat_plot <- rbind(a,b,c)
tpr2 <- ggplot(data = dat_plot, mapping = aes(x=sample.size, y=tpr, fill=sample.size), position = "dodge") +
    geom_boxplot(width = 0.5) + 
    scale_fill_manual(values = c("gray100", "gray60", "gray30")) +
    xlab("sample size, N") +
    ylab("true postive rate") + 
    ggtitle("Experiment 2 \n Meta-Feature Selection True Positive Rate") + ThemeMain; tpr2

a <- data.frame(fpr = base$fpr, sample.size = "100")
b <- data.frame(fpr = mSample$fpr, sample.size = "200")
c <- data.frame(fpr = lSample$fpr, sample.size = "500")
dat_plot <- rbind(a,b,c)
fpr2 <- ggplot(data = dat_plot, mapping = aes(x=sample.size, y=fpr, fill=sample.size), position = "dodge") +
    geom_boxplot(width = 0.5) + 
    scale_fill_manual(values = c("gray100", "gray60", "gray30")) +
    xlab("sample size, N") +
    ylab("false postive rate") + 
    ggtitle("Experiment 2 \n Meta-Feature Selection False Positive Rate") + ThemeMain; fpr2

a <- data.frame(acc = base$acc, sample.size = "100")
b <- data.frame(acc = mSample$acc, sample.size = "200")
c <- data.frame(acc = lSample$acc, sample.size = "500")
dat_plot <- rbind(a,b,c)
acc2 <- ggplot(data = dat_plot, mapping = aes(x=sample.size, y=acc, fill=sample.size), position = "dodge") +
    geom_boxplot(width = 0.5) + 
    scale_fill_manual(values = c("gray100", "gray60", "gray30")) +
    xlab("sample size, N") +
    ylab("accuracy") + 
    ggtitle("Experiment 2 \n Meta-Feature Selection Accuracy") + ThemeMain; acc2
#########################################################################################################################
# 3: large feature size vs. small, ridge-lasso, high signal
set.seed(1)
mP <- simCoxExt(N = 100, p = 500, q = 50, true_cindex = 0.8, snr_z = 2, nSim = 100)
lP <- simCoxExt(N = 100, p = 1000, q = 50, true_cindex = 0.8, snr_z = 2, nSim = 100)
a <- data.frame(cindex = base$cidx_test, feature.size="200", Method="Ridge_Lasso")
b <- data.frame(cindex = base$cidx_test_noext, feature.size="200", Method="Ridge")
c <- data.frame(cindex = mP$cidx_test, feature.size="500", Method="Ridge_Lasso")
d <- data.frame(cindex = mP$cidx_test_noext, feature.size="500", Method="Ridge")
e <- data.frame(cindex = lP$cidx_test, feature.size="1000", Method="Ridge_Lasso")
f <- data.frame(cindex = lP$cidx_test_noext, feature.size="1000", Method="Ridge")
dat_plot <- rbind(a,b,c,d,e,f)
dat_plot$feature.size <- factor(dat_plot$feature.size, levels = c("200","500","1000"))
cidx3 <- ggplot(data = dat_plot, mapping = aes(x = feature.size, y = cindex, fill = Method), position = "dodge") +
    geom_boxplot(width = 0.6) + 
    geom_hline(aes(yintercept = 0.8), linetype = "dashed") +
    scale_fill_manual(values = c("gray30", "gray80")) + 
    xlab("feature size, p") +
    ylab("test C-index") 
cidx3 <- cidx3 + ggtitle("Experiment 3 \n Prediction Performance") + ThemeMain; cidx3

a <- data.frame(tpr = base$tpr, feature.size = "200")
b <- data.frame(tpr = mP$tpr, feature.size = "500")
c <- data.frame(tpr = lP$tpr, feature.size = "1000")
dat_plot <- rbind(a,b,c)
dat_plot$feature.size <- factor(dat_plot$feature.size, levels = c("200","500","1000"))
tpr3 <- ggplot(data = dat_plot, mapping = aes(x=feature.size, y=tpr, fill=feature.size), position = "dodge") +
    geom_boxplot(width = 0.5) + 
    scale_fill_manual(values = c("gray100", "gray60", "gray30")) +
    xlab("feature size, p") +
    ylab("true postive rate") + 
    ggtitle("Experiment 3 \n Meta-Feature Selection True Positive Rate") + ThemeMain; tpr3

a <- data.frame(fpr = base$fpr, feature.size = "200")
b <- data.frame(fpr = mP$fpr, feature.size = "500")
c <- data.frame(fpr = lP$fpr, feature.size = "1000")
dat_plot <- rbind(a,b,c)
dat_plot$feature.size <- factor(dat_plot$feature.size, levels = c("200","500","1000"))
fpr3 <- ggplot(data = dat_plot, mapping = aes(x=feature.size, y=fpr, fill=feature.size), position = "dodge") +
    geom_boxplot(width = 0.5) + 
    scale_fill_manual(values = c("gray100", "gray60", "gray30")) +
    xlab("feature size, p") +
    ylab("false postive rate") +
    ggtitle("Experiment 3 \n Meta-Feature Selection False Positive Rate") + ThemeMain; fpr3

a <- data.frame(acc = base$acc, feature.size = "200")
b <- data.frame(acc = mP$acc, feature.size = "500")
c <- data.frame(acc = lP$acc, feature.size = "1000")
dat_plot <- rbind(a,b,c)
dat_plot$feature.size <- factor(dat_plot$feature.size, levels = c("200","500","1000"))
acc3 <- ggplot(data = dat_plot, mapping = aes(x=feature.size, y=acc, fill=feature.size), position = "dodge") +
    geom_boxplot(width = 0.5) + 
    scale_fill_manual(values = c("gray100", "gray60", "gray30")) +
    xlab("feature size, p") +
    ylab("accuracy") + 
    ggtitle("Experiment 3 \n Meta-Feature Selection Accuracy") + ThemeMain; acc3
#########################################################################################################################
# 4: large number of external variable vs. small, ridge-lasso, high signal
set.seed(1)
sQ <- simCoxExt(N = 100, p = 200, q = 20, n_nonzero_a = 10, true_cindex = 0.8, snr_z = 2, nSim = 100)
lQ <- simCoxExt(N = 100, p = 200, q = 100, n_nonzero_a = 10, true_cindex = 0.8, snr_z = 2, nSim = 100)
a <- data.frame(cindex = sQ$cidx_test, external.size="20", Method="Ridge_Lasso")
b <- data.frame(cindex = sQ$cidx_test_noext, external.size="20", Method="Ridge")
c <- data.frame(cindex = base$cidx_test, external.size="50", Method="Ridge_Lasso")
d <- data.frame(cindex = base$cidx_test_noext, external.size="50", Method="Ridge")
e <- data.frame(cindex = lQ$cidx_test, external.size="100", Method="Ridge_Lasso")
f <- data.frame(cindex = lQ$cidx_test_noext, external.size="100", Method="Ridge")
dat_plot <- rbind(a,b,c,d,e,f)
dat_plot$external.size <- factor(dat_plot$external.size, levels = c("20","50","100"))
cidx4 <- ggplot(data = dat_plot, mapping = aes(x = external.size, y = cindex, fill = Method), position = "dodge") +
    geom_boxplot(width = 0.6) + 
    geom_hline(aes(yintercept = 0.8), linetype = "dashed") +
    scale_fill_manual(values = c("gray30", "gray80")) +   
    xlab("meta-feature size, q") +
    ylab("test C-index") 
cidx4 <- cidx4 + ggtitle("Experiment 4 \n Prediction Performance") + ThemeMain; cidx4

a <- data.frame(tpr = sQ$tpr, external.size = "20")
b <- data.frame(tpr = base$tpr, external.size = "50")
c <- data.frame(tpr = lQ$tpr, external.size = "100")
dat_plot <- rbind(a,b,c)
dat_plot$external.size <- factor(dat_plot$external.size, levels = c("20","50","100"))
tpr4 <- ggplot(data = dat_plot, mapping = aes(x=external.size, y=tpr, fill=external.size), position = "dodge") +
    geom_boxplot(width = 0.5) + 
    scale_fill_manual(values = c("gray100", "gray60", "gray30")) +
    xlab("meta-feature size, q") +
    ylab("true postive rate") +
    ggtitle("Experiment 4 \n Meta-Feature Selection True Positive Rate") + ThemeMain; tpr4

a <- data.frame(fpr = sQ$fpr, external.size = "20")
b <- data.frame(fpr = base$fpr, external.size = "50")
c <- data.frame(fpr = lQ$fpr, external.size = "100")
dat_plot <- rbind(a,b,c)
dat_plot$external.size <- factor(dat_plot$external.size, levels = c("20","50","100"))
fpr4 <- ggplot(data = dat_plot, mapping = aes(x=external.size, y=fpr, fill=external.size), position = "dodge") +
    geom_boxplot(width = 0.5) + 
    scale_fill_manual(values = c("gray100", "gray60", "gray30")) +
    xlab("meta-feature size, q") +
    ylab("false postive rate") +
    ggtitle("Experiment 4 \n Meta-Feature Selection False Positive Rate") + ThemeMain; fpr4

a <- data.frame(acc = sQ$acc, external.size = "20")
b <- data.frame(acc = base$acc, external.size = "50")
c <- data.frame(acc = lQ$acc, external.size = "100")
dat_plot <- rbind(a,b,c)
dat_plot$external.size <- factor(dat_plot$external.size, levels = c("20","50","100"))
acc4 <- ggplot(data = dat_plot, mapping = aes(x=external.size, y=acc, fill=external.size), position = "dodge") +
    geom_boxplot(width = 0.5) + 
    scale_fill_manual(values = c("gray100", "gray60", "gray30")) +
    xlab("meta-feature size, q") +
    ylab("accuracy") + 
    ggtitle("Experiment 4 \n Meta-Feature Selection Accuracy") + ThemeMain; acc4
#########################################################################################################################
# 5: high data true C-index vs. low, ridge-lasso, high signal
set.seed(1)
sC <- simCoxExt(N = 100, p = 200, q = 50, true_cindex = 0.75, snr_z = 2, nSim = 100)
lC <- simCoxExt(N = 100, p = 200, q = 50, true_cindex = 0.85, snr_z = 2, nSim = 100)
a <- data.frame(cindex = sC$cidx_test, true.Cindex="0.75", Method="Ridge_Lasso")
b <- data.frame(cindex = sC$cidx_test_noext, true.Cindex="0.75", Method="Ridge")
c <- data.frame(cindex = base$cidx_test, true.Cindex="0.80", Method="Ridge_Lasso")
d <- data.frame(cindex = base$cidx_test_noext, true.Cindex="0.80", Method="Ridge")
e <- data.frame(cindex = lC$cidx_test, true.Cindex="0.85", Method="Ridge_Lasso")
f <- data.frame(cindex = lC$cidx_test_noext, true.Cindex="0.85", Method="Ridge")
dat_plot <- rbind(a,b,c,d,e,f)
cidx5 <- ggplot(data = dat_plot, mapping = aes(x = true.Cindex, y = cindex, fill = Method), position = "dodge") +
    geom_boxplot(width = 0.6) + 
    scale_fill_manual(values = c("gray30", "gray80")) +   
    xlab("data true C-index") +
    ylab("test C-index") 
cidx5 <- cidx5 + ggtitle("Experiment 5 \n Prediction Performance") + ThemeMain; cidx5

a <- data.frame(tpr = sC$tpr, true.Cindex = "0.75")
b <- data.frame(tpr = base$tpr, true.Cindex = "0.80")
c <- data.frame(tpr = lQ$tpr, true.Cindex = "0.85")
dat_plot <- rbind(a,b,c)
tpr5 <- ggplot(data = dat_plot, mapping = aes(x=true.Cindex, y=tpr, fill=true.Cindex), position = "dodge") +
    geom_boxplot(width = 0.5) + 
    scale_fill_manual(values = c("gray100", "gray60", "gray30")) +
    xlab("data true C-index") +
    ylab("true postive rate") +
    ggtitle("Experiment 5 \n Meta-Feature Selection True Positive Rate") + ThemeMain; tpr5

a <- data.frame(fpr = sC$fpr, true.Cindex = "0.75")
b <- data.frame(fpr = base$fpr, true.Cindex = "0.80")
c <- data.frame(fpr = lC$fpr, true.Cindex = "0.85")
dat_plot <- rbind(a,b,c)
fpr5 <- ggplot(data = dat_plot, mapping = aes(x=true.Cindex, y=fpr, fill=true.Cindex), position = "dodge") +
    geom_boxplot(width = 0.5) + 
    scale_fill_manual(values = c("gray100", "gray60", "gray30")) +
    xlab("data true C-index") +
    ylab("false postive rate") +
    ggtitle("Experiment 5 \n Meta-Feature Selection False Positive Rate") + ThemeMain; fpr5

a <- data.frame(acc = sC$acc, true.Cindex = "0.75")
b <- data.frame(acc = base$acc, true.Cindex = "0.80")
c <- data.frame(acc = lC$acc, true.Cindex = "0.85")
dat_plot <- rbind(a,b,c)
acc5 <- ggplot(data = dat_plot, mapping = aes(x=true.Cindex, y=acc, fill=true.Cindex), position = "dodge") +
    geom_boxplot(width = 0.5) + 
    scale_fill_manual(values = c("gray100", "gray60", "gray30")) +
    xlab("data true C-index") +
    ylab("accuracy") + 
    ggtitle("Experiment 5 \n Meta-Feature Selection Accuracy") + ThemeMain; acc5
#########################################################################################################################
# 6: large correlation of x vs small correlation, ridge-lasso, high signal
set.seed(1)
srho <- simCoxExt(N = 100, p = 200, q = 50, true_cindex = 0.8, corr_x = 0.2, snr_z = 2, nSim = 100)
lrho <- simCoxExt(N = 100, p = 200, q = 50, true_cindex = 0.8, corr_x = 0.8, snr_z = 2, nSim = 100)
a <- data.frame(cindex = srho$cidx_test, correlation="0.2", Method="Ridge_Lasso")
b <- data.frame(cindex = srho$cidx_test_noext, correlation="0.2", Method="Ridge")
c <- data.frame(cindex = base$cidx_test, correlation="0.5", Method="Ridge_Lasso")
d <- data.frame(cindex = base$cidx_test_noext, correlation="0.5", Method="Ridge")
e <- data.frame(cindex = lrho$cidx_test, correlation="0.8", Method="Ridge_Lasso")
f <- data.frame(cindex = lrho$cidx_test_noext, correlation="0.8", Method="Ridge")
dat_plot <- rbind(a,b,c,d,e,f)
cidx6 <- ggplot(data = dat_plot, mapping = aes(x = correlation, y = cindex, fill = Method), position = "dodge") +
    geom_boxplot(width = 0.6) + 
    geom_hline(aes(yintercept = 0.8), linetype = "dashed") +
    scale_fill_manual(values = c("gray30", "gray80")) +   
    xlab("correlation, rho") +
    ylab("test C-index") 
cidx6 <- cidx6 + ggtitle("Experiment 6 \n Prediction Performance") + ThemeMain; cidx6

a <- data.frame(tpr = srho$tpr, correlation = "0.2")
b <- data.frame(tpr = base$tpr, correlation = "0.5")
c <- data.frame(tpr = lrho$tpr, correlation = "0.8")
dat_plot <- rbind(a,b,c)
tpr6 <- ggplot(data = dat_plot, mapping = aes(x=correlation, y=tpr, fill=correlation), position = "dodge") +
    geom_boxplot(width = 0.5) + 
    scale_fill_manual(values = c("gray100", "gray60", "gray30")) +
    xlab("correlation, rho") +
    ylab("true postive rate") +
    ggtitle("Experiment 6 \n Meta-Feature Selection True Positive Rate") + ThemeMain; tpr6

a <- data.frame(fpr = srho$fpr, correlation = "0.2")
b <- data.frame(fpr = base$fpr, correlation = "0.5")
c <- data.frame(fpr = lrho$fpr, correlation = "0.8")
dat_plot <- rbind(a,b,c)
fpr6 <- ggplot(data = dat_plot, mapping = aes(x=correlation, y=fpr, fill=correlation), position = "dodge") +
    geom_boxplot(width = 0.5) + 
    scale_fill_manual(values = c("gray100", "gray60", "gray30")) +
    xlab("correlation, rho") +
    ylab("false postive rate") +
    ggtitle("Experiment 6 \n Meta-Feature Selection False Positive Rate") + ThemeMain; fpr6

a <- data.frame(acc = srho$acc, correlation = "0.2")
b <- data.frame(acc = base$acc, correlation = "0.5")
c <- data.frame(acc = lrho$acc, correlation = "0.8")
dat_plot <- rbind(a,b,c)
acc6 <- ggplot(data = dat_plot, mapping = aes(x=correlation, y=acc, fill=correlation), position = "dodge") +
    geom_boxplot(width = 0.5) + 
    scale_fill_manual(values = c("gray100", "gray60", "gray30")) +
    xlab("correlation, rho") +
    ylab("accuracy") + 
    ggtitle("Experiment 6 \n Meta-Feature Selection Accuracy") + ThemeMain; acc6
#########################################################################################################################
save(cidx1, cidx2, cidx3, cidx4, cidx5, cidx6, tpr1, tpr2, tpr3, tpr4, tpr5, tpr6, fpr1, fpr2, fpr3, fpr4, fpr5, fpr6,
     acc1, acc2, acc3, acc4, acc5, acc6,
     file = "simCox_plot.RData")
save(zeroext, base, medext, lowext, mSample, lSample, mP, lP, sQ, lQ, sC, lC, srho, lrho, 
     file = "simCox.RData")

load("simCox_plot.RData")

library(gridExtra)
pdf("simCox_plot.pdf", width = 16, height = 9)
grid.arrange(cidx1, cidx2, ncol=2)
grid.arrange(cidx3, cidx4, ncol=2)
grid.arrange(cidx5, cidx6, ncol=2)
grid.arrange(acc1, acc2, ncol=2)
grid.arrange(acc3, acc4, ncol=2)
grid.arrange(acc5, acc6, ncol=2)
grid.arrange(tpr1, tpr2, ncol=2)
grid.arrange(tpr3, tpr4, ncol=2)
grid.arrange(tpr5, tpr6, ncol=2)
grid.arrange(fpr1, fpr2, ncol=2)
grid.arrange(fpr3, fpr4, ncol=2)
grid.arrange(fpr5, fpr6, ncol=2)
dev.off()

cidx1 <- cidx1 + ggtitle("(a). Experiment 1 \n Prediction Performance") + ThemeMain
cidx2 <- cidx2 + ggtitle("(b). Experiment 2 \n Prediction Performance") + ThemeMain
cidx3 <- cidx3 + ggtitle("(c). Experiment 3 \n Prediction Performance") + ThemeMain
cidx4 <- cidx4 + ggtitle("(d). Experiment 4 \n Prediction Performance") + ThemeMain
png("sim1.png", units="in", width=18, height=12, res=300)
grid.arrange(cidx1, cidx2, cidx3, cidx4, ncol=2)
dev.off()

tpr1 = tpr1 + ggtitle("(a). Meta-Feature Selection \n True Positive Rate")
fpr1 = fpr1 + ggtitle("(b). Meta-Feature Selection \n False Positive Rate")
png("sim2.png", units="in", width=18, height=8, res=300)
grid.arrange(tpr1, fpr1, ncol=2)
dev.off()




