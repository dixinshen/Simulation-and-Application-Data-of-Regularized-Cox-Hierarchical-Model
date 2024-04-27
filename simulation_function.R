# simulate function
simCoxExt <- function(N, p, q, true_cindex, snr_z, nSim, corr_x = 0.5, cov_x = "auto", alpha1 = 0, alpha2 = 1, 
                      nonzero_a = 0.2, n_nonzero_a = NULL, alpha_same_sign = TRUE, effect_size = 0.2, external = NULL) {
    shape <- 8
    scale <- 5
    # fix alphas
    a <- rep(0, q)
    if (!is.null(n_nonzero_a)) {
        a[1:n_nonzero_a] <- effect_size
    } else {
        if (alpha_same_sign == TRUE) {
            a[c(1:floor(nonzero_a*q))] <- effect_size
        } else {
            a[c(1:floor(nonzero_a*q))] <- c(rep(effect_size,floor(nonzero_a*q/2)), rep(-effect_size, (floor(nonzero_a*q)-floor(nonzero_a*q/2))))
        }
    }
    
    if ( is.null(external) ) {
        z <- rbinom(p*q, size = 1, prob = 0.1)
        z <- matrix(z, nrow = p, ncol = q)
        mode(z) <- "double"
    } else {
        z <- external
    }
    
    # simulate betas
    if (alpha1==0) {
        var_za <- var(drop(z %*% a))  # empirical variance
        var_b <- var_za / snr_z
        b <- drop(z %*% a) + sqrt(var_b)*rnorm(p)
    } else if (alpha1==1) {
        var_za <- var(drop(z %*% a))  # empirical variance
        var_b <- var_za / snr_z
        b <- drop(z %*% a) + extraDistr::rlaplace(p, 0, sqrt(var_b/2))
    }
    
    
    # set parameters
    if (cov_x == "auto") {
        sigma <- matrix(NA, nrow = p, ncol = p)
        for (j in 1:p) {
            for (i in 1:p) {
                sigma[i, j] <- corr_x^abs(i - j)
            }
        }
    } else if (cov_x == "cs") {
        block <- matrix( rep(.5, (p/q)^2), nrow=p/q )
        diag(block) <- 1
        sigma <- as.matrix(bdiag(replicate(q, block, simplify=F)))
    }
    x <- mvrnorm((1e4+N), rep(0, p), sigma)
    
    # fix true c index
    for (stdev in seq(0.01, 10, 0.01)) {
        t <- scale * (-log(runif(1e4+N)) * exp(-drop(x%*%b)))^(1/shape) + rnorm((1e4+N), sd = stdev)
        t[t<=0] <- min(t[t>0])
        t[t>20] <- 20    # follow up time is 20
        c <- rexp((1e4+N), 0.06)     # produce ratio of event to censor at around 3/1
        c[c>20] <- 20
        time <- pmin(t,c)
        status <- ifelse(t<c, 1, 0)
        y <- cbind(time, status)
        if ( abs(Cindex(drop(x%*%b), y) - true_cindex) < 0.001 ) break
    }
    cat("Data truc C-index is ", Cindex(drop(x%*%b), y), "\n", sep = "")
    x_val <- x[1:N, ]
    y_val <- y[1:N, ]
    x_test <- x[-(1:N), ]
    y_test <- y[-(1:N), ]
    
    cidx_test_r_l <- numeric()
    coefs_r_l <- mat.or.vec(nSim, p)
    cidx_test_l_l <- numeric()
    coefs_l_l <- mat.or.vec(nSim, p)
    cidx_test_ridge <- numeric()
    coefs_ridge <- mat.or.vec(nSim, p)
    cidx_test_lasso <- numeric()
    coefs_lasso <- mat.or.vec(nSim, p)
    tpr <- numeric()
    fpr <- numeric()
    acc <- numeric()
    npos <- numeric()
    alphas_r_l <- mat.or.vec(nSim, q)
    alphas_l_l <- mat.or.vec(nSim, q)
    for (j in 1:nSim) {
        # xrnet
        x <- mvrnorm(N, rep(0, p), sigma)
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
        
        fit_r_l <- xrnet(x, y, z, family="cox", 
                     penalty_main=define_ridge(penalty_ratio=0.01), 
                     penalty_external=define_lasso(penalty_ratio=0.01),
                     control=list(max_iterations=1e6))

        fit_l_l <- xrnet(x, y, z, family="cox", 
                     penalty_main=define_lasso(penalty_ratio=0.01), 
                     penalty_external=define_lasso(penalty_ratio=0.01),
                     control=list(max_iterations=1e6))
        
        beta_r_l <- fit_r_l$betas
        alpha_r_l <- fit_r_l$alphas
        cidx_val_r_l <- apply(beta_r_l, c(2,3), function(l) Cindex(drop(x_val%*%l), y_val))
        wh_r_l <- drop(which(cidx_val_r_l==max(cidx_val_r_l), arr.ind = T))
        coefs_r_l[j, ] <- beta_r_l[, wh_r_l[1], wh_r_l[2]]
        alphas_r_l[j, ] <- alpha_r_l[, wh_r_l[1], wh_r_l[2]]
        tpr[j] <- sum(alphas_r_l[j, 1:floor(nonzero_a*q)]!=0) / floor(nonzero_a*q)
        fpr[j] <- sum(alphas_r_l[j, -(1:floor(nonzero_a*q))]!=0) / (q - floor(nonzero_a*q))
        acc[j] <- sum(ifelse(a==0, 0, 1) == ifelse(alphas_r_l[j, ]==0, 0, 1)) / q
        npos[j] <- sum(alphas_r_l[j, ] != 0)
        cidx_test_r_l[j] <- Cindex(drop(x_test%*%beta_r_l[, wh_r_l[1], wh_r_l[2]]), y_test)
        
        beta_l_l <- fit_l_l$betas
        alpha_l_l <- fit_l_l$alphas
        cidx_val_l_l <- apply(beta_l_l, c(2,3), function(l) Cindex(drop(x_val%*%l), y_val))
        wh_l_l <- drop(which(cidx_val_l_l==max(cidx_val_l_l), arr.ind = T))
        coefs_l_l[j, ] <- beta_l_l[, wh_l_l[1], wh_l_l[2]]
        alphas_l_l[j, ] <- alpha_l_l[, wh_l_l[1], wh_l_l[2]]
        cidx_test_l_l[j] <- Cindex(drop(x_test%*%beta_l_l[, wh_l_l[1], wh_l_l[2]]), y_test)
        
        # ridge
        fit_ridge <- glmnet(x, y, alpha = 0, family = "cox")
        beta_ridge <- as.matrix(fit_ridge$beta)
        cidx_val_ridge <- apply(beta_ridge, 2, function(l) Cindex(drop(x_val%*%l), y_val))
        wh_ridge <- which.max(cidx_val_ridge)
        coefs_ridge[j, ] <- beta_ridge[, wh_ridge]
        cidx_test_ridge[j] <- Cindex(drop(x_test%*%beta_ridge[, wh_ridge]), y_test)
        
        # lasso
        fit_lasso <- glmnet(x, y, alpha = 1, family = "cox")
        beta_lasso <- as.matrix(fit_lasso$beta)
        cidx_val_lasso <- apply(beta_lasso, 2, function(l) Cindex(drop(x_val%*%l), y_val))
        wh_lasso <- which.max(cidx_val_lasso)
        coefs_lasso[j, ] <- beta_ridge[, wh_lasso]
        cidx_test_lasso[j] <- Cindex(drop(x_test%*%beta_lasso[, wh_lasso]), y_test)
    }
    
    return(list(cidx_test_r_l=cidx_test_r_l, cidx_test_l_l=cidx_test_l_l,
                cidx_test_ridge=cidx_test_ridge, cidx_test_lasso=cidx_test_lasso,
                betas_r_l=coefs_r_l, betas_l_l=coefs_l_l, 
                betas_ridge=coefs_ridge, betas_lasso=coefs_lasso,
                alphas=alphas_r_l, tpr=tpr, fpr=fpr, acc=acc, number_positive=npos))
}
