# Regularized Cox function, linearized computation of weights and working responses

cdl_Cox <- function(y, x, z = NULL, standardize = TRUE, external = TRUE, alpha = 1, alpha1 = 0, alpha2 = 1, thresh = 1e-7,
                    iter_max = 100) {
    ordered <- order(y[,1])
    y <- y[ordered, ]
    x <- as.matrix(x[ordered, ])
    
    n <- NROW(y)
    p <- NCOL(x)
    w <- rep(1/n, n)
    
    # standardize x
    if (standardize == TRUE) { 
        xm <- colMeans(x)
        x_norm <- sweep(x, 2, xm, "-")
        xs <- drop(sqrt(crossprod(w, x_norm^2)))
        x_norm <- sweep(x_norm, 2, xs, "/")
    } else {
        x_norm <- x
    }
    
    if (external == FALSE) {
        nlam <- 100
        # initialize beta, weights, working response
        b_prime <- rep(0, p)
        betaHats <- mat.or.vec(nr=nlam, nc=p)
        
        t_event <- y[,1][y[,2]==1]   # event time
        D <- unique(t_event)   # unique event time
        m <- length(D)   # number of unique times
        d <- numeric()
        for(i in 1:m) {
            d[i] <- sum(t_event==D[i])
        }
        
        Ck_prime <- 0
        Ck <- numeric()
        Ri <- numeric()
        for (k in 1:n) {
            Ck[k] <- Ck_prime
            for (j in (Ck_prime+1):m) {
                if (D[j] <= y[k,1]) {
                    Ck[k] <- Ck[k] + 1
                    Ri[Ck[k]] <- n - k + 1
                } else {
                    break
                }
                Ck_prime <- Ck[k]
            }
            if (Ck_prime==m) {break}
        }
        if (k < n) {Ck[(k+1):n] <- Ck_prime}
        Ri <- c(Ri, 0)
        Ck <- c(0, Ck)
        
        eta <- as.vector(x_norm %*% b_prime)
        exp_eta <- exp(eta)
        sum_exp_eta_prime <- 0
        sum_exp_eta <- numeric()
        for (i in m:1) {
            sum_exp_eta[i] <- sum_exp_eta_prime + sum( exp_eta[(n-Ri[i]+1):(n-Ri[i+1])] )
            sum_exp_eta_prime <- sum_exp_eta[i]
        }
        u_prime <- 0
        u2_prime <- 0
        u <- numeric()
        u2 <- numeric()
        for (k in 1:n) {
            if (Ck[k+1] - Ck[k] == 0) {
                u[k] <- u_prime 
                u2[k] <- u2_prime
            } else {
                u[k] <- u_prime + d[Ck[k+1]] / sum_exp_eta[Ck[k+1]]
                u2[k] <- u2_prime + d[Ck[k+1]] / (sum_exp_eta[Ck[k+1]])^2
            }
            u_prime <- u[k]
            u2_prime <- u2[k]
        }
        W <- exp_eta*u - exp_eta*exp_eta*u2   # weights
        r <- w * (W*eta + y[,2] - exp_eta*u)  # residuals = weights * working response / n
        
        # compute penalty path
        if (alpha > 0) {
            lambdaMax <- max( abs(crossprod(x_norm, r)) ) / alpha
        } else if (alpha == 0) {
            lambdaMax <- 1000 * max( abs(crossprod(x_norm, r)) )
        }
        if (n >= p) {
            lambdaMin <- 0.0001*lambdaMax
            logLambda <- seq(log(lambdaMax), log(lambdaMin), length.out = nlam)
            lambdas <- exp(logLambda)
        } else {
            lambdaMin <- 0.01*lambdaMax
            logLambda <- seq(log(lambdaMax), log(lambdaMin), length.out = nlam)
            lambdas <- exp(logLambda)
        }
        
        b_current <- b_prime
        # penalty path loop
        for(l in 1:nlam) {
            lambda_current <- lambdas[l] 
            dev_out <- 1e10
            iter <- 0
            
            # outer reweighted least square loop
            while(dev_out >= thresh) {
                b_old <- b_current 
                x2w <- drop(crossprod(x_norm^2, W*w))
                
                # inner coordinate descent loop
                converge_inner <- FALSE
                while (!converge_inner) {
                    iter <- iter + 1
                    dev_in <- 0.0
                    for(j in 1:p) {
                        bj <- b_current[j]
                        rj <- r + w*W*(b_current[j]*x_norm[,j])
                        wls <- sum(x_norm[,j]*rj)
                        arg <- abs(wls) - alpha*lambda_current
                        if (arg > 0.0) {
                            b_current[j] <- sign(wls)*arg / (x2w[j] + (1-alpha)*lambda_current)
                        } else {
                            b_current[j] <- 0
                        }
                        del <- b_current[j] - bj
                        if (abs(del) > 0.0) {
                            dev_in <- max(dev_in, abs(del))
                            r <- r - del*(x_norm[,j])*W*w
                        }
                    }
                    if (dev_in < thresh) {converge_inner <- TRUE}
                }
                dev_out <- max( abs(b_current-b_old) ) 
                
                # update weights and working responses
                eta <- x_norm %*% b_current
                exp_eta <- exp(eta)
                sum_exp_eta_prime <- 0
                sum_exp_eta <- numeric()
                for (i in m:1) {
                    sum_exp_eta[i] <- sum_exp_eta_prime + sum( exp_eta[(n-Ri[i]+1):(n-Ri[i+1])] )
                    sum_exp_eta_prime <- sum_exp_eta[i]
                }
                u_prime <- 0
                u2_prime <- 0
                u <- numeric()
                u2 <- numeric()
                for (k in 1:n) {
                    if (Ck[k+1] - Ck[k] == 0) {
                        u[k] <- u_prime 
                        u2[k] <- u2_prime
                    } else {
                        u[k] <- u_prime + d[Ck[k+1]] / sum_exp_eta[Ck[k+1]]
                        u2[k] <- u2_prime + d[Ck[k+1]] / (sum_exp_eta[Ck[k+1]])^2
                    }
                    u_prime <- u[k]
                    u2_prime <- u2[k]
                }
                W <- exp_eta*u - exp_eta*exp_eta*u2  
                r <- w * (W*eta + y[,2] - exp_eta*u) - w*W*eta 
            } # outer while loop
            
            if (standardize == TRUE) {
                betaHats[l, ] <- b_current/xs
            } else {
                betaHats[l, ] <- b_current
            }
        }
        return(list(coef=betaHats, lambda=lambdas))
    }
    
    # Incorporating external information
    if (external == TRUE) {
        q <- NCOL(z)
        # standardize xz
        if (standardize == TRUE) {
            zm <- colMeans(z)
            z_norm <- sweep(z, 2, zm, "-")
            zs <- drop(sqrt(crossprod(rep(1/p, p), z_norm^2)))
            z_norm <- sweep(z, 2, zs, "/")
            xz <- x_norm %*% z_norm
            xzs <- drop(sqrt(crossprod(w, sweep(x%*%z, 2, colMeans(x%*%z), "-")^2)))
            xzs_norm <- drop(sqrt(crossprod(w, sweep(xz, 2, colMeans(xz), "-")^2)))
        } else {
            xz <- x %*% z
        }
        
        g_l11 <- rep(0, p)
        a_l11 <- rep(0, q)
        nlam <- 20
        betaHats <- array(0, dim = c(p+q, nlam, nlam))
        
        t_event <- y[,1][y[,2]==1]   # event time
        D <- unique(t_event)   # unique event time
        m <- length(D)   # number of unique times
        d <- numeric()
        for(i in 1:m) {
            d[i] <- sum(t_event==D[i])
        }
        
        Ck_prime <- 0
        Ck <- numeric()
        Ri <- numeric()
        for (k in 1:n) {
            Ck[k] <- Ck_prime
            for (j in (Ck_prime+1):m) {
                if (D[j] <= y[k,1]) {
                    Ck[k] <- Ck[k] + 1
                    Ri[Ck[k]] <- n - k + 1
                } else {
                    break
                }
                Ck_prime <- Ck[k]
            }
            if (Ck_prime==m) {break}
        }
        if (k < n) {Ck[(k+1):n] <- Ck_prime}
        Ri <- c(Ri, 0)
        Ck <- c(0, Ck)
        
        eta <- as.vector(x_norm %*% g_l11 + xz %*% a_l11)
        exp_eta <- exp(eta)
        sum_exp_eta_prime <- 0
        sum_exp_eta <- numeric()
        for (i in m:1) {
            sum_exp_eta[i] <- sum_exp_eta_prime + sum( exp_eta[(n-Ri[i]+1):(n-Ri[i+1])] )
            sum_exp_eta_prime <- sum_exp_eta[i]
        }
        u_prime <- 0
        u2_prime <- 0
        u <- numeric()
        u2 <- numeric()
        for (k in 1:n) {
            if (Ck[k+1] - Ck[k] == 0) {
                u[k] <- u_prime 
                u2[k] <- u2_prime
            } else {
                u[k] <- u_prime + d[Ck[k+1]] / sum_exp_eta[Ck[k+1]]
                u2[k] <- u2_prime + d[Ck[k+1]] / (sum_exp_eta[Ck[k+1]])^2
            }
            u_prime <- u[k]
            u2_prime <- u2[k]
        }
        W_l11 <- exp_eta*u - exp_eta*exp_eta*u2   # weights
        r_l11 <- w * (W_l11*eta + y[,2] - exp_eta*u)  # residuals = weights * working response / n
        
        # compute penalty pathm 
        if (alpha2 > 0) {
            lambda2_max <- max( abs(crossprod(xz, r_l11)) ) / alpha2
        } else if (alpha2 == 0) {
            lambda2_max <- 1000 * max( abs(crossprod(xz, r_l11)) ) 
        }
        if (alpha1 > 0 ) {
            lambda1_max <- max( abs(crossprod(x_norm, r_l11)) ) / alpha1
        } else if (alpha1 == 0) {
            lambda1_max <- 1000 * max( abs(crossprod(x_norm, r_l11)) )
        }
        if (n >= (p+q)) {
            lambda2_min <- 0.0001*lambda2_max
            lambda1_min <- 0.0001*lambda1_max
            lambda2 <- exp( seq(log(lambda2_max), log(lambda2_min), length.out = nlam) )
            lambda1 <- exp( seq(log(lambda1_max), log(lambda1_min), length.out = nlam) )
        } else {
            lambda2_min <- 0.01*lambda2_max
            lambda1_min <- 0.01*lambda1_max
            lambda2 <- exp( seq(log(lambda2_max), log(lambda2_min), length.out = nlam) )
            lambda1 <- exp( seq(log(lambda1_max), log(lambda1_min), length.out = nlam) )
        }
        
        # penalty path loop
        for (l2 in 1:nlam) {
            lambda2_current <- lambda2[l2]
            g_current <- g_l11
            a_current <- a_l11
            W <- W_l11
            r <- r_l11
            
            for (l1 in 1:nlam) {
                lambda1_current <- lambda1[l1]
                dev_out <- 1e10
                iter_outer <- 0
                
                # outer reweighted least square loop
                while (dev_out >= thresh) {
                    g_old <- g_current
                    a_old <- a_current
                    x2w <- drop(crossprod(x_norm^2, W*w))
                    xz2w <- drop(crossprod(xz^2, W*w))
                    
                    # inner coordinate descent loop 
                    iter_inner <- 0
                    converge_inner <- FALSE
                    while (!converge_inner) {
                        criteria0 <- sum(r^2/(2*w*W)) + sum(g_current^2)*lambda1_current/2 + sum(abs(a_current))*lambda2_current
                        dev_in <- 0.0
                        for (j in 1:p) {
                            gj <- g_current[j]
                            rj <- r + w*W*(g_current[j]*x_norm[,j])
                            wls <- sum(x_norm[,j]*rj)
                            arg <- abs(wls) - alpha1*lambda1_current
                            if (arg>0.0) {
                                g_current[j] <- sign(wls)*arg / (x2w[j] + (1-alpha1)*lambda1_current)
                            } else {
                                g_current[j] <- 0
                            }
                            del <- g_current[j] - gj
                            if(abs(del) > 0.0) {
                                dev_in <- max(dev_in, abs(del))
                                r <- r - del*(x_norm[,j])*W*w
                            }
                        }
                        for (k in 1:q) {
                            ak <- a_current[k]
                            rk <- r + w*W*(a_current[k]*xz[,k])
                            wls <- sum(xz[,k]*rk)
                            arg <- abs(wls) - alpha2*lambda2_current
                            if(arg>0.0) {
                                a_current[k] <- sign(wls)*arg / (xz2w[k] + (1-alpha2)*lambda2_current)
                            } else {
                                a_current[k] <- 0
                            }
                            del <- a_current[k] - ak
                            if(abs(del) > 0.0) {
                                dev_in <- max(dev_in, abs(del))
                                r <- r - del*(xz[,k])*W*w
                            }
                        }
                        if (dev_in < thresh) {converge_inner <- TRUE}
                        iter_inner <- iter_inner + 1
                        if (dev_in > 10) {
                            # cat("estimation blow out, break out from inner coordinate descent loop\n")
                            break
                        }
                        if (iter_inner > iter_max) {
                            # cat("Inner coordinate descent loop max iteration reached, dev_in: ", dev_in, "\n")
                            break
                        }
                    } # inner while loop
                    
                    dev_out <- max( abs(c(g_current,a_current) - c(g_old,a_old)) )
                    iter_outer <- iter_outer + 1
                    if(dev_out > 10) {
                        # cat("estimation blow out, break out from outer reweighted loop\n")
                        break
                    }
                    if (iter_outer > iter_max) {
                        # cat("Outer reweighted loop max iteration reached, dev_out: ", dev_out, "\n")
                        break
                    }
                    
                    # update weights and working responses
                    eta <- x_norm %*% g_current + xz %*% a_current
                    exp_eta <- exp(eta)
                    sum_exp_eta_prime <- 0
                    sum_exp_eta <- numeric()
                    for (i in m:1) {
                        sum_exp_eta[i] <- sum_exp_eta_prime + sum( exp_eta[(n-Ri[i]+1):(n-Ri[i+1])] )
                        sum_exp_eta_prime <- sum_exp_eta[i]
                    }
                    u_prime <- 0
                    u2_prime <- 0
                    u <- numeric()
                    u2 <- numeric()
                    for (k in 1:n) {
                        if (Ck[k+1] - Ck[k] == 0) {
                            u[k] <- u_prime 
                            u2[k] <- u2_prime
                        } else {
                            u[k] <- u_prime + d[Ck[k+1]] / sum_exp_eta[Ck[k+1]]
                            u2[k] <- u2_prime + d[Ck[k+1]] / (sum_exp_eta[Ck[k+1]])^2
                        }
                        u_prime <- u[k]
                        u2_prime <- u2[k]
                    }
                    W <- exp_eta*u - exp_eta*exp_eta*u2   # weights
                    r <- w * (W*eta + y[,2] - exp_eta*u) - w*W*eta  # residuals = weights * working response / n
                } # outer while loop
                
                if (dev_in > 10 | dev_out > 10) {
                    cat("break out from inner lambda1 loop\n")
                    break
                }
                
                # keep the residual of lambda1 at l1=1 for warm start
                if (l1==1) {
                    g_l11 <- g_current
                    a_l11 <- a_current
                    W_l11 <- W
                    r_l11 <- r
                }
                
                if (standardize == TRUE) {
                    b <- g_current + drop(z_norm %*% a_current)
                    betaHats[,l1,l2] <- c(b, a_current) / c(xs, xzs/xzs_norm) 
                } else {
                    b <- g_current + drop(z %*% a_current)
                    betaHats[,l1,l2] <- c(b, a_current)
                }
                
            } # inner for lambda1 loop
            
            # create estimate names, column names
            est_names <- c(paste("beta_", 1:p, sep=""), paste("alpha_", 1:q, sep=""))
            dimnames(betaHats)[[1]] <- est_names
            
            if (dev_in > 10 | dev_out > 10) {
                cat("break out from outer lambda2 loop\n")
                break
            }
        } # outer for lambda2 loop
        
        return(list(coef=betaHats, lambda=rbind(lambda1, lambda2)))
    }
    
}

