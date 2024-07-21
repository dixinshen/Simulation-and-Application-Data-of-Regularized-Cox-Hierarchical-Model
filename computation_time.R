shape <- 8
scale <- 5
N = 100
q = 50
# fix alphas
a <- rep(0, q)
a[1:10] <- 0.2

for (p in c(200, 1000, 10000)) {
    z <- rbinom(p*q, size = 1, prob = 0.1)
    z <- matrix(z, nrow = p, ncol = q)
    mode(z) <- "double"
    
    var_za <- var(drop(z %*% a))  # empirical variance
    var_b <- var_za / 2
    b <- drop(z %*% a) + sqrt(var_b)*rnorm(p)
    
    sigma <- matrix(NA, nrow = p, ncol = p)
    for (j in 1:p) {
        for (i in 1:p) {
            sigma[i, j] <- 0.5^abs(i - j)
        }
    }
    x <- mvnfast::rmvn(1000, rep(0, p), sigma)
    # fix true c index
    for (stdev in seq(0.01, 10, 0.01)) {
        t <- scale * (-log(runif(1000)) * exp(-drop(x%*%b)))^(1/shape) + rnorm(1000, sd = stdev)
        t[t<=0] <- min(t[t>0])
        t[t>20] <- 20    # follow up time is 20
        c <- rexp(1000, 0.06)     # produce ratio of event to censor at around 3/1
        c[c>20] <- 20
        time <- pmin(t,c)
        status <- ifelse(t<c, 1, 0)
        y <- cbind(time, status)
        if ( abs(Cindex(drop(x%*%b), y) - 0.8) < 0.001 ) break
    }
    
    x <- mvnfast::rmvn(N, rep(0, p), sigma)
    t <- scale * (-log(runif(N)) * exp(-drop(x%*%b)))^(1/shape) + rnorm(N, sd = stdev)
    t[t<=0] <- min(t[t>0])
    t[t>20] <- 20    # follow up time is 20
    c <- rexp(N, 0.06)     # produce ratio of event to censor at around 7/3
    c[c>20] <- 20
    time <- pmin(t,c)
    status <- ifelse(t<c, 1, 0)
    table(status)
    y <- cbind(time, status)
    
    ordered <- order(y[,1])
    y <- y[ordered, ]
    x <- as.matrix(x[ordered, ])
    
    start <- Sys.time()
    for (j in 1:100) {
        fit_r_l <- xrnet(x, y, z, family="cox", 
                         penalty_main=define_ridge(penalty_ratio=0.01), 
                         penalty_external=define_lasso(penalty_ratio=0.01),
                         control=list(max_iterations=1e6))
    }
    end <- Sys.time()
    cat("The average computational time of N=100, p=", p, ", q=50, is ", (end-start)/100, "\n")
}

start <- Sys.time()
fit_r_l <- xrnet(x, y, z, family="cox", 
                 penalty_main=define_ridge(penalty_ratio=0.01), 
                 penalty_external=define_lasso(penalty_ratio=0.01),
                 control=list(max_iterations=1e6))
end <- Sys.time()
end - start














