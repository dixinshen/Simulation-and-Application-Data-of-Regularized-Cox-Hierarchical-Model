library(MASS)
library(glmnet)
library(xrnet)
library(ggplot2)
library(ggpubr)

source("simulation_function.R")
#####################################################################################

#####################################################################################
ThemeMain <- theme(
    plot.title = element_text(size=12, face='bold', hjust=0.5),
)

# base case: N = 100, p = 200, q = 50, true_cindex = 0.8, snr_z = 2, cov_x = "auto", alpha1 = 0, alpha2 = 1
# 1: high external signal vs. low, ridge-lasso
set.seed(2024)
start <- Sys.time()
zeroext <- simCoxExt(N = 100, p = 200, q = 50, true_cindex = 0.8, snr_z = 5e-2, nSim = 100)
base <- simCoxExt(N = 100, p = 200, q = 50, true_cindex = 0.8, snr_z = 2, nSim = 100)
medext <- simCoxExt(N = 100, p = 200, q = 50, true_cindex = 0.8, snr_z = 1, nSim = 100)
lowext <- simCoxExt(N = 100, p = 200, q = 50, true_cindex = 0.8, snr_z = 0.1, nSim = 100)
end <- Sys.time()
end-start
a <- data.frame(cindex = base$cidx_test_r_l, signal=2, Method="Ridge_Lasso")
b <- data.frame(cindex = base$cidx_test_l_l, signal=2, Method="Lasso_Lasso")
c <- data.frame(cindex = base$cidx_test_ridge, signal=2, Method="Ridge")
d <- data.frame(cindex = base$cidx_test_lasso, signal=2, Method="Lasso")
e <- data.frame(cindex = medext$cidx_test_r_l, signal = 1, Method="Ridge_Lasso")
f <- data.frame(cindex = medext$cidx_test_l_l, signal = 1, Method="Lasso_Lasso")
g <- data.frame(cindex = medext$cidx_test_ridge, signal=1, Method="Ridge")
h <- data.frame(cindex = medext$cidx_test_lasso, signal=1, Method="Lasso")
i <- data.frame(cindex = lowext$cidx_test_r_l, signal=0.1, Method="Ridge_Lasso")
j <- data.frame(cindex = lowext$cidx_test_l_l, signal=0.1, Method="Lasso_Lasso")
k <- data.frame(cindex = lowext$cidx_test_ridge, signal=0.1, Method="Ridge")
l <- data.frame(cindex = lowext$cidx_test_lasso, signal=0.1, Method="Lasso")
m <- data.frame(cindex = zeroext$cidx_test_r_l, signal=0, Method="Ridge_Lasso")
n <- data.frame(cindex = zeroext$cidx_test_l_l, signal=0, Method="Lasso_Lasso")
o <- data.frame(cindex = zeroext$cidx_test_ridge, signal=0, Method="Ridge")
r <- data.frame(cindex = zeroext$cidx_test_lasso, signal=0, Method="Lasso")
dat_plot <- rbind(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,r)
dat_plot$signal <- as.factor(dat_plot$signal)

cidx11 <- ggplot(data = dat_plot[!dat_plot$Method%in%c("Lasso","Lasso_Lasso"),], 
                 mapping = aes(x = signal, y = cindex, fill = Method), 
                 position = "dodge") +
    geom_boxplot(width = 0.5) + 
    geom_hline(aes(yintercept = 0.8), linetype = "dashed") +
    scale_fill_manual(values = c("gray30", "gray80")) +
    xlab("Meta-feature signal-noise ratio, SNR") +
    ylab("Test C-index")
cidx11 <- cidx11 + ggtitle("Experiment 1") + ThemeMain; cidx11

cidx12 <- ggplot(data = dat_plot, mapping = aes(x = signal, y = cindex, fill = Method), position = "dodge") +
    geom_boxplot(width = 0.5) + 
    geom_hline(aes(yintercept = 0.8), linetype = "dashed") +
    scale_fill_manual(values = c("gray10", "gray70", "gray40", "gray100")) +
    xlab("Meta-feature signal-noise ratio, SNR") +
    ylab("Test C-index")
cidx12 <- cidx12 + ggtitle("Experiment 1") + ThemeMain; cidx12

a <- data.frame(tpr = base$tpr, signal = 2)
b <- data.frame(tpr = medext$tpr, signal = 1)
c <- data.frame(tpr = lowext$tpr, signal = 0.1)
d <- data.frame(tpr = zeroext$tpr, signal = 0)
dat_plot <- rbind(a,b,c,d)
dat_plot$signal <- as.factor(dat_plot$signal)
tpr1 <- ggplot(data = dat_plot, mapping = aes(x=signal, y=tpr, fill=signal), position = "dodge") +
    geom_boxplot(width = 0.5) + 
    scale_fill_manual(values = c("gray100", "gray70", "gray40", "gray10")) +
    xlab("Meta-feature signal-noise ratio, SNR") +
    ylab("True postive rate") + 
    ggtitle("Meta-feature selection") + ThemeMain; tpr1

a <- data.frame(fpr = base$fpr, signal = 2)
b <- data.frame(fpr = medext$fpr, signal = 1)
c <- data.frame(fpr = lowext$fpr, signal = 0.1)
d <- data.frame(fpr = zeroext$fpr, signal = 0)
dat_plot <- rbind(a,b,c,d)
dat_plot$signal <- as.factor(dat_plot$signal)
fpr1 <- ggplot(data = dat_plot, mapping = aes(x=signal, y=fpr, fill=signal), position = "dodge") +
    geom_boxplot(width = 0.5) + 
    scale_fill_manual(values = c("gray100", "gray70", "gray40", "gray10")) +
    xlab("Meta-feature signal-noise ratio, SNR") +
    ylab("False postive rate") + 
    ggtitle("Meta-feature selection") + ThemeMain; fpr1
#########################################################################################################################
# 2: large sample size vs. small, ridge-lasso, high signal
set.seed(2024)
mSample <- simCoxExt(N = 200, p = 200, q = 50, true_cindex = 0.8, snr_z = 2, nSim = 100)
lSample <- simCoxExt(N = 500, p = 200, q = 50, true_cindex = 0.8, snr_z = 2, nSim = 100)
a <- data.frame(cindex = base$cidx_test_r_l, sample.size="100", Method="Ridge_Lasso")
b <- data.frame(cindex = base$cidx_test_l_l, sample.size="100", Method="Lasso_Lasso")
c <- data.frame(cindex = base$cidx_test_ridge, sample.size="100", Method="Ridge")
d <- data.frame(cindex = base$cidx_test_lasso, sample.size="100", Method="Lasso")
e <- data.frame(cindex = mSample$cidx_test_r_l, sample.size="200", Method="Ridge_Lasso")
f <- data.frame(cindex = mSample$cidx_test_l_l, sample.size="200", Method="Lasso_Lasso")
g <- data.frame(cindex = mSample$cidx_test_ridge, sample.size="200", Method="Ridge")
h <- data.frame(cindex = mSample$cidx_test_lasso, sample.size="200", Method="Lasso")
i <- data.frame(cindex = lSample$cidx_test_r_l, sample.size="500", Method="Ridge_Lasso")
j <- data.frame(cindex = lSample$cidx_test_l_l, sample.size="500", Method="Lasso_Lasso")
k <- data.frame(cindex = lSample$cidx_test_ridge, sample.size="500", Method="Ridge")
l <- data.frame(cindex = lSample$cidx_test_lasso, sample.size="500", Method="Lasso")
dat_plot <- rbind(a,b,c,d,e,f,g,h,i,j,k,l)
cidx21 <- ggplot(data = dat_plot[!dat_plot$Method%in%c("Lasso","Lasso_Lasso"),], 
                 mapping = aes(x = sample.size, y = cindex, fill = Method), 
                 position = "dodge") +
    geom_boxplot(width = 0.6) + 
    geom_hline(aes(yintercept = 0.8), linetype = "dashed") +
    scale_fill_manual(values = c("gray30", "gray80")) +    
    xlab("Sample size, N") +
    ylab("Test C-index") 
cidx21 <- cidx21 + ggtitle("Experiment 2") + ThemeMain; cidx21

cidx22 <- ggplot(data = dat_plot, mapping = aes(x = sample.size, y = cindex, fill = Method), position = "dodge") +
    geom_boxplot(width = 0.5) + 
    geom_hline(aes(yintercept = 0.8), linetype = "dashed") +
    scale_fill_manual(values = c("gray10", "gray70", "gray40", "gray100")) +
    xlab("Sample size, N") +
    ylab("Test C-index")
cidx22 <- cidx22 + ggtitle("Experiment 2") + ThemeMain; cidx22
#########################################################################################################################
# 3: large feature size vs. small, ridge-lasso, high signal
set.seed(2024)
mP <- simCoxExt(N = 100, p = 500, q = 50, true_cindex = 0.8, snr_z = 2, nSim = 100)
lP <- simCoxExt(N = 100, p = 1000, q = 50, true_cindex = 0.8, snr_z = 2, nSim = 100)
a <- data.frame(cindex = base$cidx_test_r_l, feature.size="200", Method="Ridge_Lasso")
b <- data.frame(cindex = base$cidx_test_l_l, feature.size="200", Method="Lasso_Lasso")
c <- data.frame(cindex = base$cidx_test_ridge, feature.size="200", Method="Ridge")
d <- data.frame(cindex = base$cidx_test_lasso, feature.size="200", Method="Lasso")
e <- data.frame(cindex = mP$cidx_test_r_l, feature.size="500", Method="Ridge_Lasso")
f <- data.frame(cindex = mP$cidx_test_l_l, feature.size="500", Method="Lasso_Lasso")
g <- data.frame(cindex = mP$cidx_test_ridge, feature.size="500", Method="Ridge")
h <- data.frame(cindex = mP$cidx_test_lasso, feature.size="500", Method="Lasso")
i <- data.frame(cindex = lP$cidx_test_r_l, feature.size="1000", Method="Ridge_Lasso")
j <- data.frame(cindex = lP$cidx_test_l_l, feature.size="1000", Method="Lasso_Lasso")
k <- data.frame(cindex = lP$cidx_test_ridge, feature.size="1000", Method="Ridge")
l <- data.frame(cindex = lP$cidx_test_lasso, feature.size="1000", Method="Lasso")
dat_plot <- rbind(a,b,c,d,e,f,g,h,i,j,k,l)
dat_plot$feature.size <- factor(dat_plot$feature.size, levels = c("200","500","1000"))
cidx31 <- ggplot(data = dat_plot[!dat_plot$Method%in%c("Lasso","Lasso_Lasso"),], 
                 mapping = aes(x = feature.size, y = cindex, fill = Method), 
                 position = "dodge") +
    geom_boxplot(width = 0.6) + 
    geom_hline(aes(yintercept = 0.8), linetype = "dashed") +
    scale_fill_manual(values = c("gray30", "gray80")) +    
    xlab("Feature size, p") +
    ylab("Test C-index") 
cidx31 <- cidx31 + ggtitle("Experiment 3") + ThemeMain; cidx31

cidx32 <- ggplot(data = dat_plot, mapping = aes(x = feature.size, y = cindex, fill = Method), position = "dodge") +
    geom_boxplot(width = 0.5) + 
    geom_hline(aes(yintercept = 0.8), linetype = "dashed") +
    scale_fill_manual(values = c("gray10", "gray70", "gray40", "gray100")) +
    xlab("Feature size, p") +
    ylab("Test C-index")
cidx32 <- cidx32 + ggtitle("Experiment 3") + ThemeMain; cidx32
#########################################################################################################################
# 4: large number of external variable vs. small, ridge-lasso, high signal
set.seed(2024)
sQ <- simCoxExt(N = 100, p = 200, q = 20, true_cindex = 0.8, snr_z = 2, nSim = 100)
lQ <- simCoxExt(N = 100, p = 200, q = 100, true_cindex = 0.8, snr_z = 2, nSim = 100)
a <- data.frame(cindex = sQ$cidx_test_r_l, external.size="20", Method="Ridge_Lasso")
b <- data.frame(cindex = sQ$cidx_test_l_l, external.size="20", Method="Lasso_Lasso")
c <- data.frame(cindex = sQ$cidx_test_ridge, external.size="20", Method="Ridge")
d <- data.frame(cindex = sQ$cidx_test_lasso, external.size="20", Method="Lasso")
e <- data.frame(cindex = base$cidx_test_r_l, external.size="50", Method="Ridge_Lasso")
f <- data.frame(cindex = base$cidx_test_l_l, external.size="50", Method="Lasso_Lasso")
g <- data.frame(cindex = base$cidx_test_ridge, external.size="50", Method="Ridge")
h <- data.frame(cindex = base$cidx_test_lasso, external.size="50", Method="Lasso")
i <- data.frame(cindex = lQ$cidx_test_r_l, external.size="100", Method="Ridge_Lasso")
j <- data.frame(cindex = lQ$cidx_test_l_l, external.size="100", Method="Lasso_Lasso")
k <- data.frame(cindex = lQ$cidx_test_ridge, external.size="100", Method="Ridge")
l <- data.frame(cindex = lQ$cidx_test_lasso, external.size="100", Method="Lasso")
dat_plot <- rbind(a,b,c,d,e,f,g,h,i,j,k,l)
dat_plot$external.size <- factor(dat_plot$external.size, levels = c("20","50","100"))
cidx41 <- ggplot(data = dat_plot[!dat_plot$Method%in%c("Lasso","Lasso_Lasso"),], 
                 mapping = aes(x = external.size, y = cindex, fill = Method), 
                 position = "dodge") +
    geom_boxplot(width = 0.6) + 
    geom_hline(aes(yintercept = 0.8), linetype = "dashed") +
    scale_fill_manual(values = c("gray30", "gray80")) +    
    xlab("Meta-feature size, q") +
    ylab("Test C-index") 
cidx41 <- cidx41 + ggtitle("Experiment 4") + ThemeMain; cidx41

cidx42 <- ggplot(data = dat_plot, mapping = aes(x = external.size, y = cindex, fill = Method), position = "dodge") +
    geom_boxplot(width = 0.5) + 
    geom_hline(aes(yintercept = 0.8), linetype = "dashed") +
    scale_fill_manual(values = c("gray10", "gray70", "gray40", "gray100")) +
    xlab("Meta-feature size, p") +
    ylab("Test C-index")
cidx42 <- cidx42 + ggtitle("Experiment 4") + ThemeMain; cidx42
#########################################################################################################################

#########################################################################################################################
save(cidx11, cidx12, cidx21, cidx22, cidx31, cidx32, cidx41, cidx42, tpr1, fpr1,
     file = "simCox_plot.RData")
save(zeroext, base, medext, lowext, mSample, lSample, mP, lP, sQ, lQ,
     file = "simCox.RData")

load("simCox_plot.RData")

g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}

mylegend <- g_legend(cidx11)

cidx11 <- cidx11 + theme(legend.position="none")
cidx21 <- cidx21 + theme(legend.position="none")
cidx31 <- cidx31 + theme(legend.position="none")
cidx41 <- cidx41 + theme(legend.position="none")

sim1 <- ggarrange(cidx11,cidx21,cidx31,cidx41, 
                  nrow=2, ncol=2, common.legend = TRUE, legend="bottom")
ggsave("sim1.pdf", units="mm", width=300, height=200, dpi=300)

tpr1 = tpr1 + theme(legend.position="none")
fpr1 = fpr1 + theme(legend.position="none")

sim2 <- ggarrange(tpr1, fpr1, nrow=1, ncol=2, common.legend=T, legend="bottom")
ggsave("sim2.pdf", units="mm", width=300, height=130, dpi=300)

sim3 <- ggarrange(cidx12,cidx22,cidx32,cidx42, 
                  nrow=2, ncol=2, common.legend = TRUE, legend="bottom")
ggsave("sim3.pdf", units="mm", width=300, height=200, dpi=300)
#########################################################################################################################

#########################################################################################################################




