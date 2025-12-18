setwd("/home/hoanguc3m/Dropbox/WP11/Code/htvpAU/")
#install.packages("fattvpVAR_0.2.0.tar.gz", repos = NULL, type = "source")

library(readxl)
library(Matrix)
library(fattvpVAR)
library(profvis)
library(invgamma)

#Data_AUS <- read_excel("Data AUS.xlsx", sheet = "Q")
Data_AUS <- read_excel("Data251207.xlsx", sheet = "Q")

vars::VARselect(y = Data_AUS[2:115,2:4], lag.max = 10)

K <- 3 # 3 series
p <- 2 # number of lags

dataraw <- Data_AUS[2:115,2:4]
dataraw[,1] <- dataraw[,1] /100

y <- data.matrix(dataraw[(p+1):nrow(dataraw),])
y0 <- data.matrix(dataraw[1:p,])

n.ahead = 20
{
  library(ggplot2)
  library(gridExtra)
  library(ggthemr)
  ggthemr('fresh')
  get_irf_tvp_stdS <- function(Chain, impulse.variable = 1, response.variable = 2,
                               atT = NULL, n.ahead = 24, draw.plot = FALSE){
    # TODO: Check IRF in case of fat tail.
    K <- ncol(Chain$data$y)
    p <- Chain$data$p

    k_alp <- K*(K-1)/2 # dimension of the impact matrix
    k_beta <- K^2*p + K # number of VAR coefficients
    k_beta_div_K <- K*p + 1
    is_tv <- Chain$data$inits$is_tv
    n_tv <- sum(is_tv)# # of time-varying equations

    k_a_eq <- seq(0, K-1)
    id_a <- cbind(cumsum(k_a_eq) - k_a_eq + 1,cumsum(k_a_eq))
    count_seqa <- list(); count_seqb <- list()
    count_seqa[[1]] <- 0; count_seqb[[1]] <- seq(1, k_beta_div_K)
    for (ii in c(2:K)){
      count_seqa[[ii]] <- seq(id_a[ii,1], id_a[ii,2])
      count_seqb[[ii]] <- ((ii-1)*k_beta_div_K+1):(ii*k_beta_div_K)
    }

    idx_b_tv <- (kronecker(is_tv,matrix(1,nrow = K*p+1,ncol = 1))==1)   # index for time-varying betas
    idx_a_tv <- matrix(FALSE, nrow = k_alp, ncol = 1)            # construct index for time-varying alphas
    for (ii in 2:K){
      if (is_tv[ii] == 1){
        idx_a_tv[ count_seqa[[ii]]] <- TRUE
      }
    }


    ndraws <- nrow(Chain$store_alp0)
    t_max <- nrow(Chain$data$y)

    if (is.null(atT)){
      atT = t_max
    }


    dist <- Chain$data$dist
    SV <- TRUE #Chain$data$priors


    if (SV) {
      sig <- exp(Chain$store_h[, atT, ]*0.5)
    }

    if (dist == "Student"){
      W <- sqrt(Chain$store_w[, atT,])
    }


    H.chol <- array(diag(K), dim = c(K,K, n.ahead+1))
    beta <- get_post(Chain,element = "beta")
    alp <- get_post(Chain,element = "alpha")

    beta <- beta[,,-c( seq(from = 1, to = k_beta, by = k_beta_div_K) )] # Remove constant


    out <- matrix(NA, ndraws, n.ahead + 1)
    out_all <- array(NA, c(K,K,ndraws, n.ahead + 1))

    # Check stable draws
    stable_idx <- rep(FALSE, ndraws)
    for (j in c(1:ndraws)){
      A0 <- matrix(0, nrow = K, ncol = K)
      A0[lower.tri(A0)] <- alp[j,atT,]
      diag(A0) <- 1
      A_inv <- solve(A0)
      B_tmp = matrix(beta[j,atT,], nrow = K, byrow = T)
      Ainv_B <- B_tmp
      for (i in c(1:p)){
        Ainv_B[,(K* (i-1) + 1):(K*i)] <- A_inv %*% B_tmp[,(K* (i-1) + 1):(K*i)]
      }

      if (stab_mat(B = Ainv_B, K = K, p = p) ){
        stable_idx[j] <- TRUE
      }
    }
    draw_stab <- sample((1:ndraws)[stable_idx], size = ndraws, replace = TRUE)

    # Compute IR for each MC iteration
    for (j in c(1:ndraws)){
      samp <- draw_stab[j]
      # Cholesky of VCV matrix, depending on specification

      A0 <- matrix(0, nrow = K, ncol = K)
      A0[lower.tri(A0)] <- alp[samp,atT,]
      diag(A0) <- 1
      A_inv <- solve(A0)

      s = n.ahead+1
      sigtemp <- sig[samp,c(1:K)]
      if (dist == "Gaussian"){
        H.chol[,,s] <- A_inv %*% diag(sigtemp)
      }
      if (dist == "Student" ){
        #wtemp <- W[samp,c(1:K)]
        # H.chol[,,s] <- A_inv %*% diag(wtemp) %*% diag(sigtemp)
        H.chol[,,s] <- A_inv %*% diag(sigtemp)
      }
      for (s in c(1:n.ahead)){
        H.chol[,,s] <- H.chol[,,n.ahead+1]
      }

      B_tmp = matrix(beta[samp,atT,], nrow = K, byrow = T)
      Ainv_B <- B_tmp
      for (i in c(1:p)){
        Ainv_B[,(K* (i-1) + 1):(K*i)] <- A_inv %*% B_tmp[,(K* (i-1) + 1):(K*i)]
      }

      # Compute Impulse Responses
      # if (stab_mat(B = Ainv_B, K = K, p = p) ){
      aux <- IRFmats(B = Ainv_B,
                     H.chol = H.chol)

      out[j,] <- aux[response.variable, impulse.variable, ]
      out_all[,,j,] <- aux
      # }

    }

    return(list(irf_all = out,
                mean_irf = apply(out, MARGIN = 2, FUN = mean),
                mean_all = apply(out_all, MARGIN = c(1,2,4), FUN = mean) )  )
  }

  impres2dplotv2 <- function(tmp, name = "", cid = "", yrange = c(0,1)){
    mean_irf <- apply(tmp$irf_all, MARGIN = 2, FUN = mean)
    q10_irf <- apply(tmp$irf_all, MARGIN = 2, FUN = quantile, probs = 0.1)
    q90_irf <- apply(tmp$irf_all, MARGIN = 2, FUN = quantile, probs = 0.9)

    xrange <- (0:(length(mean_irf)-1))

    data_nu <- data.frame(Time = xrange,
                          irf = mean_irf)

    data_nu1 <- data.frame(Time = c(xrange, rev(xrange)) ,
                           irfCI = c(q10_irf, rev(q90_irf)))

    p1 <- ggplot() + geom_line(data=data_nu, mapping=aes(x=Time, y=irf), color = "#ff0c00") +
      geom_polygon(data=data_nu1, mapping=aes(x=Time, y=irfCI), fill = "#5ba6d6", alpha = 0.5) +
      ggtitle(name) + xlab("") + ylab(cid) + theme(plot.title = element_text(hjust = 0.5, size=11)) # + ylim(yrange)
    return(p1)
  }

  plot_3x3 <- function(Chain = G000_obj){
    atT <- nrow(Chain$data$y)
    p1 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = Chain,
                                                impulse.variable = 1, response.variable = 1,
                                                n.ahead = n.ahead, atT = atT),
                         name = "EPU", cid = "EPU") ; p1

    p2 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = Chain,
                                                impulse.variable = 1, response.variable = 2,
                                                n.ahead = n.ahead, atT = atT),
                         name = "", cid = "Unemployment") ; p2
    p3 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = Chain,
                                                impulse.variable = 1, response.variable = 3,
                                                n.ahead = n.ahead, atT = atT),
                         name = "", cid = "GDP growth") ; p3

    p4 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = Chain,
                                                impulse.variable = 2, response.variable = 1,
                                                n.ahead = n.ahead, atT = atT),
                         name = "Unemployment") ; p4


    p5 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = Chain,
                                                impulse.variable = 2, response.variable = 2,
                                                n.ahead = n.ahead, atT = atT),
                         name = "") ; p5

    p6 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = Chain,
                                                impulse.variable = 2, response.variable = 3,
                                                n.ahead = n.ahead, atT = atT),
                         name = "") ; p6

    p7 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = Chain,
                                                impulse.variable = 3, response.variable = 1,
                                                n.ahead = n.ahead, atT = atT),
                         name = "GDP growth") ; p7

    p8 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = Chain,
                                                impulse.variable = 3, response.variable = 2,
                                                n.ahead = n.ahead, atT = atT),
                         name = "") ; p8
    p9 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = Chain,
                                                impulse.variable = 3, response.variable = 3,
                                                n.ahead = n.ahead, atT = atT),
                         name = "") ; p9

    #png(filename='/home/hoanguc3m/Dropbox/WP11/Code/htvpAU/img/IR/IRF_G000_std_2025Q1.png', width = 1200, height = 1200)

    #par(mar=c(2,1,3,5))
    return(grid.arrange( p1, p2, p3,
                         p4, p5, p6,
                         p7, p8, p9,
                         nrow = 3, ncol = 3, as.table = FALSE))
    #dev.off()
  }
  plot_2x2 <- function(Chain = G000_obj, nameID){
    atT <- nrow(Chain$data$y)
    p1 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = Chain,
                                                impulse.variable = 1, response.variable = 1,
                                                n.ahead = n.ahead, atT = atT),
                         name = nameID[1], cid = nameID[1]) ; p1

    p2 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = Chain,
                                                impulse.variable = 1, response.variable = 2,
                                                n.ahead = n.ahead, atT = atT),
                         name = "", cid = nameID[2]) ; p2

    p4 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = Chain,
                                                impulse.variable = 2, response.variable = 1,
                                                n.ahead = n.ahead, atT = atT),
                         name = nameID[2]) ; p4


    p5 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = Chain,
                                                impulse.variable = 2, response.variable = 2,
                                                n.ahead = n.ahead, atT = atT),
                         name = "") ; p5


    return(grid.arrange( p1, p2,
                         p4, p5,
                         nrow = 2, ncol = 2, as.table = FALSE))

  }
}

# Full data
priors <- get_prior_minnesota(y = y, p = p, intercept=TRUE)
inits <- list(samples = 10000, burnin = 2000, thin = 2)
inits$is_tv = c(0,0,0); G000_obj <- fitTVPGaussSV_ASIS(y, y0, p, priors, inits)
Plot1_G000SVFull <- plot_3x3(Chain = G000_obj)
ggsave("img/IR/IRF_Fig1.png", Plot1_G000SVFull)

# Short data
priors <- get_prior_minnesota(y = y[1:90,], p = p, intercept=TRUE)
inits <- list(samples = 10000, burnin = 2000, thin = 2)
inits$is_tv = c(0,0,0); G000Short_obj <- fitTVPGaussSV_ASIS(y[1:90,], y0, p, priors, inits)

Plot1_G000SVShort <- plot_3x3(Chain = G000Short_obj)
ggsave("img/IR/IRF_Fig2_Short.png", Plot1_G000SVShort)

# Full data with tight prior
priors <- get_prior_minnesota(y = y, p = p, intercept=TRUE)
priors$V_a0_prior <- c(0.1,0.1,0.1)
inits <- list(samples = 10000, burnin = 2000, thin = 2)
inits$is_tv = c(0,0,0); G000tight_obj <- fitTVPGaussSV_ASIS(y, y0, p, priors, inits)
Plot1_G000SVTight <- plot_3x3(Chain = G000tight_obj)
ggsave("img/IR/IRF_Fig3_tight.png", Plot1_G000SVTight)

# Full data with G000 non SV
priors <- get_prior_minnesota(y = y, p = p, intercept=TRUE)
inits <- list(samples = 10000, burnin = 2000, thin = 2)
inits$is_tv = c(0,0,0); G000nonSV_obj <- fitTVPGaussnonSV_ASIS(y, y0, p, priors, inits)
Plot1_G000nonSV <- plot_3x3(Chain = G000nonSV_obj)
ggsave("img/IR/IRF_Fig4_nonSV.png", Plot1_G000nonSV)


# Full data with G000 non SV
priors <- get_prior_minnesota(y = y, p = p, intercept=TRUE, lambda2 = 0.05)
priors$hyper_ab <- list(hyper_a = rep(0.0001, 3),
                        hyper_b = rep(0.0001,21))
inits <- list(samples = 10000, burnin = 2000, thin = 2)
inits$is_tv = c(1,1,1); G111nonSV_obj <- fitTVPGaussnonSV_ASIS(y, y0, p, priors, inits)
Plot1_G111nonSV <- plot_3x3(Chain = G111nonSV_obj)
ggsave("img/IR/IRF_Fig10TVPatight_nonSV.png", Plot1_G111nonSV)
plot(colMeans(G111nonSV_obj$store_alp[,,3]) )
plot(colMeans(G111nonSV_obj$store_beta[,,4]) )
hist(G111nonSV_obj$store_alp0[,3])
inits$is_tv = c(1,1,1); G111_SV_obj <- fitTVPGaussSV_ASIS(y, y0, p, priors, inits)
Plot1_G000nonSV <- plot_3x3(Chain = G111_SV_obj)

# Bivar data
Serie_name <- c("EPU", "Unemployment","GDP Growth")

priors <- get_prior_minnesota(y = y[,c(1,2)], p = p, intercept=TRUE)
inits <- list(samples = 10000, burnin = 2000, thin = 2)
inits$is_tv = c(0,0); G000SV12_obj <- fitTVPGaussSV_ASIS(y[,c(1,2)], y0[,c(1,2)], p, priors, inits)
Plot2_G000_obj <- plot_2x2(Chain = G000SV12_obj, nameID = Serie_name[c(1,2)])
ggsave("img/IR/IRF_Fig5_SV12.png", Plot2_G000_obj)

priors <- get_prior_minnesota(y = y[,c(1,3)], p = p, intercept=TRUE)
inits <- list(samples = 10000, burnin = 2000, thin = 2)
inits$is_tv = c(0,0); G000SV13_obj <- fitTVPGaussSV_ASIS(y[,c(1,3)], y0[,c(1,3)], p, priors, inits)
Plot2_G000_obj <- plot_2x2(Chain = G000SV13_obj, nameID = Serie_name[c(1,3)])
ggsave("img/IR/IRF_Fig6_SV13.png", Plot2_G000_obj)

priors <- get_prior_minnesota(y = y[,c(2,3)], p = p, intercept=TRUE)
inits <- list(samples = 10000, burnin = 2000, thin = 2)
inits$is_tv = c(0,0); G000SV23_obj <- fitTVPGaussSV_ASIS(y[,c(2,3)], y0[,c(2,3)], p, priors, inits)
Plot2_G000_obj <- plot_2x2(Chain = G000SV23_obj, nameID = Serie_name[c(2,3)])
ggsave("img/IR/IRF_Fig7_SV23.png", Plot2_G000_obj)

priors <- get_prior_minnesota(y = y[,c(1,3)], p = p, intercept=TRUE)
inits <- list(samples = 10000, burnin = 2000, thin = 2)
inits$is_tv = c(1,1); G000TVPSV13_obj <- fitTVPGaussSV_ASIS(y[,c(1,3)], y0[,c(1,3)], p, priors, inits)
Plot2_G000_obj <- plot_2x2(Chain = G000TVPSV13_obj, nameID = Serie_name[c(1,3)])
ggsave("img/IR/IRF_Fig8_TVPSV13.png", Plot2_G000_obj)


priors <- get_prior_minnesota(y = y[,c(1,3)], p = p, intercept=TRUE)
priors$V_a0_prior <- c(0.1,0.1,0.1)
inits <- list(samples = 10000, burnin = 2000, thin = 2)
inits$is_tv = c(1,1); G000TVPSV13t_obj <- fitTVPGaussSV_ASIS(y[,c(1,3)], y0[,c(1,3)], p, priors, inits)
Plot2_G000_obj <- plot_2x2(Chain = G000TVPSV13t_obj, nameID = Serie_name[c(1,3)])
ggsave("img/IR/IRF_Fig9_TVPSV13tight.png", Plot2_G000_obj)


# Full data with tight prior
priors <- get_prior_minnesota(y = y, p = p, intercept=TRUE, lambda2 = 0.05)
priors$V_a0_prior <- c(0.1,0.1,0.1)
inits <- list(samples = 10000, burnin = 2000, thin = 2)
inits$is_tv = c(0,0,0); G000tight_obj <- fitTVPGaussSV_ASIS(y, y0, p, priors, inits)
Plot1_G000SVTight <- plot_3x3(Chain = G000tight_obj)
ggsave("img/IR/IRF_Fig3_tighttight.png", Plot1_G000SVTight)


# Full data with loose loose prior
priors <- get_prior_minnesota(y = y, p = p, intercept=TRUE, lambda2 = 0.05)
priors$V_a0_prior <- c(100,100,100)
priors$V_b_prior <- rep(100,21)
inits <- list(samples = 10000, burnin = 2000, thin = 2)
inits$is_tv = c(0,0,0); G000tight_obj <- fitTVPGaussnonSV_ASIS(y, y0, p, priors, inits)
Plot1_G000SVTight <- plot_3x3(Chain = G000tight_obj)
ggsave("img/IR/IRF_Fig4b_looseloose.png", Plot1_G000SVTight)
