library(fattvpVAR)



{
  library(ggplot2)
  library(gridExtra)
  library(ggthemr)
  ggthemr('fresh')

  impres2dplotv2 <- function(tmp, name = "", cid = "", yrange = c(0,1)){
    mean_irf <- apply(tmp$irf_all, MARGIN = 2, FUN = mean)
    q10_irf <- apply(tmp$irf_all, MARGIN = 2, FUN = quantile, probs = 0.025)
    q90_irf <- apply(tmp$irf_all, MARGIN = 2, FUN = quantile, probs = 0.975)

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

  impres2dplotv3 <- function(tmp, name = "", cid = "", yrange = c(0,1)){
    mean_irf <- apply(tmp$irf_all, MARGIN = 2, FUN = mean)
    q10_irf <- apply(tmp$irf_all, MARGIN = 2, FUN = quantile, probs = 0.025)
    q90_irf <- apply(tmp$irf_all, MARGIN = 2, FUN = quantile, probs = 0.975)

    xrange <- (0:(length(mean_irf)-1))

    data_nu <- data.frame(Time = xrange[2:length(xrange)],
                          irf = mean_irf[2:length(xrange)])

    data_nu1 <- data.frame(Time = c(xrange, rev(xrange)) ,
                           irfCI = c(q10_irf, rev(q90_irf)))

    p1 <- ggplot() + geom_line(data=data_nu, mapping=aes(x=Time, y=irf), color = "#ff0c00") +
      geom_polygon(data=data_nu1, mapping=aes(x=Time, y=irfCI), fill = "#5ba6d6", alpha = 0.5) +
      ggtitle(name) + xlab("") + ylab(cid) + theme(plot.title = element_text(hjust = 0.5, size=11)) # + ylim(yrange)
    return(p1)
  }
}

load("/home/hoanguc3m/Downloads/htvpAU/T000.RData")

Chain <- T000_obj
atT <- nrow(T000_obj$data$y) # max
n.ahead <- 100

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

get_irf_tvp_unitS <- function(Chain, impulse.variable = 1, response.variable = 2,
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
      H.chol[,,s] <- diag(1, nrow = K, ncol = K)
    }
    if (dist == "Student" ){
      H.chol[,,s] <- diag(1, nrow = K, ncol = K)
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

RhpcBLASctl::blas_set_num_threads(1)

T000_IR_12_stdS <- parallel::mclapply(c(1:nrow(T000_obj$data$y)),
                                      FUN = function(atT) {
                                        get_irf_tvp_stdS(Chain = T000_obj,
                                                         impulse.variable = 1,
                                                         response.variable = 2,
                                                         n.ahead = n.ahead,
                                                         atT = atT)$mean_all
                                      },
                                      mc.cores = 16)

T000_IR_12_unitS <- parallel::mclapply(c(1:nrow(T000_obj$data$y)),
                                       FUN = function(atT) {
                                         get_irf_tvp_unitS(Chain = T000_obj,
                                                           impulse.variable = 1,
                                                           response.variable = 2,
                                                           n.ahead = n.ahead,
                                                           atT = atT)$mean_all
                                       },
                                       mc.cores = 16)


save(T000_IR_12_stdS,
     T000_IR_12_unitS,
     file = "/home/hoanguc3m/Downloads/htvpAU/Impulse/T000_models_IRF_new.RData")


t_max <- nrow(T000_obj$data$y)
Time <- seq(as.Date("1997/06/01"), as.Date("2025/03/01"), "quarter")
K <- 3

T000_stdS_12_irf_matrix <- array(NA, dim = c(K,K,n.ahead+1, t_max))
T000_unitS_12_irf_matrix <- array(NA, dim = c(K,K,n.ahead+1, t_max))

T_range <- 1:t_max
D_range <- tail(Time, t_max)
for (atT in c(1:t_max)){
  T000_stdS_12_irf_matrix[,,,atT] <- T000_IR_12_stdS[[atT]]
  T000_unitS_12_irf_matrix[,,,atT] <- T000_IR_12_unitS[[atT]]
}


library(plotly)
setwd("/home/hoanguc3m/Dropbox/WP11/Code/htvpAU/img/IR")

##################################

fig <- plot_ly(x = D_range,
               y = c(0:n.ahead),
               z = T000_unitS_12_irf_matrix[1,1,,T_range]) %>% add_surface() %>% layout(title = "",scene = list(xaxis = list(title = ""),yaxis = list(title = ""),zaxis = list(title = "") ))

fig <- plot_ly(x = D_range,
               y = c(0:n.ahead),
               z = T000_unitS_12_irf_matrix[1,1,,T_range]) %>% add_surface() %>% layout(title = "",scene = list(xaxis = list(title = ""),yaxis = list(title = ""),zaxis = list(title = "") ))
htmlwidgets::saveWidget(fig, file = "Model-T000-SV-UnitS-Imp1-Res1.html")

fig <- plot_ly(x = D_range,
               y = c(0:n.ahead),
               z = T000_unitS_12_irf_matrix[2,1,,T_range]) %>% add_surface() %>% layout(title = "",scene = list(xaxis = list(title = ""),yaxis = list(title = ""),zaxis = list(title = "") ))
htmlwidgets::saveWidget(fig, file = "Model-T000-SV-UnitS-Imp1-Res2.html")

fig <- plot_ly(x = D_range,
               y = c(0:n.ahead),
               z = T000_unitS_12_irf_matrix[3,1,,T_range]) %>% add_surface() %>% layout(title = "",scene = list(xaxis = list(title = ""),yaxis = list(title = ""),zaxis = list(title = "") ))
htmlwidgets::saveWidget(fig, file = "Model-T000-SV-UnitS-Imp1-Res3.html")

####


setwd("/home/hoanguc3m/Dropbox/WP11/Code/htvpAU/img/IR")

# T000_unitS_12_irf_matrix <- T000_unitS_12_irf_matrix[,,1:(1+n.ahead),]

fig <- plot_ly(x = D_range,
               y = c(0:n.ahead),
               z = T000_unitS_12_irf_matrix[1,3,,T_range]) %>% add_surface() %>% layout(title = "",scene = list(xaxis = list(title = ""),yaxis = list(title = ""),zaxis = list(title = "") ))
htmlwidgets::saveWidget(fig, file = "Model-T000-SV-UnitS-Imp3-Res1.html")

fig <- plot_ly(x = D_range,
               y = c(0:n.ahead),
               z = T000_unitS_12_irf_matrix[2,3,,T_range]) %>% add_surface() %>% layout(title = "",scene = list(xaxis = list(title = ""),yaxis = list(title = ""),zaxis = list(title = "") ))
htmlwidgets::saveWidget(fig, file = "Model-T000-SV-UnitS-Imp3-Res2.html")

fig <- plot_ly(x = D_range,
               y = c(0:n.ahead),
               z = T000_unitS_12_irf_matrix[3,3,,T_range]) %>% add_surface() %>% layout(title = "",scene = list(xaxis = list(title = ""),yaxis = list(title = ""),zaxis = list(title = "") ))
htmlwidgets::saveWidget(fig, file = "Model-T000-SV-UnitS-Imp3-Res3.html")



fig <- plot_ly(x = D_range,
               y = c(0:n.ahead),
               z = T000_unitS_12_irf_matrix[1,2,,T_range]) %>% add_surface() %>% layout(title = "",scene = list(xaxis = list(title = ""),yaxis = list(title = ""),zaxis = list(title = "") ))
htmlwidgets::saveWidget(fig, file = "Model-T000-SV-UnitS-Imp2-Res1.html")

fig <- plot_ly(x = D_range,
               y = c(0:n.ahead),
               z = T000_unitS_12_irf_matrix[2,2,,T_range]) %>% add_surface() %>% layout(title = "",scene = list(xaxis = list(title = ""),yaxis = list(title = ""),zaxis = list(title = "") ))
htmlwidgets::saveWidget(fig, file = "Model-T000-SV-UnitS-Imp2-Res2.html")

fig <- plot_ly(x = D_range,
               y = c(0:n.ahead),
               z = T000_unitS_12_irf_matrix[3,2,,T_range]) %>% add_surface() %>% layout(title = "",scene = list(xaxis = list(title = ""),yaxis = list(title = ""),zaxis = list(title = "") ))
htmlwidgets::saveWidget(fig, file = "Model-T000-SV-UnitS-Imp2-Res3.html")

###############################
{
  library(ggplot2)
  library(gridExtra)
  library(ggthemr)
  ggthemr('fresh')

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
}

#Model T000
n.ahead <- 50
atT <- 110
D_range[atT]
p1 <- impres2dplotv2(tmp = get_irf_tvp_unitS(Chain = T000_obj,
                                             impulse.variable = 1, response.variable = 1,
                                             n.ahead = n.ahead, atT = atT),
                     name = "EPU", cid = "EPU") ; p1

p2 <- impres2dplotv2(tmp = get_irf_tvp_unitS(Chain = T000_obj,
                                             impulse.variable = 1, response.variable = 2,
                                             n.ahead = n.ahead, atT = atT),
                     name = "", cid = "Unemployment") ; p2
p3 <- impres2dplotv2(tmp = get_irf_tvp_unitS(Chain = T000_obj,
                                             impulse.variable = 1, response.variable = 3,
                                             n.ahead = n.ahead, atT = atT),
                     name = "", cid = "GDP growth") ; p3

p4 <- impres2dplotv2(tmp = get_irf_tvp_unitS(Chain = T000_obj,
                                             impulse.variable = 2, response.variable = 1,
                                             n.ahead = n.ahead, atT = atT),
                     name = "Unemployment") ; p4


p5 <- impres2dplotv2(tmp = get_irf_tvp_unitS(Chain = T000_obj,
                                             impulse.variable = 2, response.variable = 2,
                                             n.ahead = n.ahead, atT = atT),
                     name = "") ; p5

p6 <- impres2dplotv2(tmp = get_irf_tvp_unitS(Chain = T000_obj,
                                             impulse.variable = 2, response.variable = 3,
                                             n.ahead = n.ahead, atT = atT),
                     name = "") ; p6

p7 <- impres2dplotv2(tmp = get_irf_tvp_unitS(Chain = T000_obj,
                                             impulse.variable = 3, response.variable = 1,
                                             n.ahead = n.ahead, atT = atT),
                     name = "GDP growth") ; p7

p8 <- impres2dplotv2(tmp = get_irf_tvp_unitS(Chain = T000_obj,
                                             impulse.variable = 3, response.variable = 2,
                                             n.ahead = n.ahead, atT = atT),
                     name = "") ; p8
p9 <- impres2dplotv2(tmp = get_irf_tvp_unitS(Chain = T000_obj,
                                             impulse.variable = 3, response.variable = 3,
                                             n.ahead = n.ahead, atT = atT),
                     name = "") ; p9

png(filename='/home/hoanguc3m/Dropbox/WP11/Code/htvpAU/img/IR/IRF_T000_unit_2025Q1.png', width = 1200*0.6, height = 1200*0.6)
#plot_grid(p12, bottom_row, labels = c('', ''), ncol = 1)
par(mar=c(2,1,3,5))
grid.arrange( p1, p2, p3,
              p4, p5, p6,
              p7, p8, p9,
              nrow = 3, ncol = 3, as.table = FALSE)
dev.off()

#######################################
#Model G000
n.ahead <- 50
atT <- 110
D_range[atT]
p1 <- impres2dplotv2(tmp = get_irf_tvp_unitS(Chain = G000_obj,
                                             impulse.variable = 1, response.variable = 1,
                                             n.ahead = n.ahead, atT = atT),
                     name = "EPU", cid = "EPU") ; p1

p2 <- impres2dplotv2(tmp = get_irf_tvp_unitS(Chain = G000_obj,
                                             impulse.variable = 1, response.variable = 2,
                                             n.ahead = n.ahead, atT = atT),
                     name = "", cid = "Unemployment") ; p2
p3 <- impres2dplotv2(tmp = get_irf_tvp_unitS(Chain = G000_obj,
                                             impulse.variable = 1, response.variable = 3,
                                             n.ahead = n.ahead, atT = atT),
                     name = "", cid = "GDP growth") ; p3

p4 <- impres2dplotv2(tmp = get_irf_tvp_unitS(Chain = G000_obj,
                                             impulse.variable = 2, response.variable = 1,
                                             n.ahead = n.ahead, atT = atT),
                     name = "Unemployment") ; p4


p5 <- impres2dplotv2(tmp = get_irf_tvp_unitS(Chain = G000_obj,
                                             impulse.variable = 2, response.variable = 2,
                                             n.ahead = n.ahead, atT = atT),
                     name = "") ; p5

p6 <- impres2dplotv2(tmp = get_irf_tvp_unitS(Chain = G000_obj,
                                             impulse.variable = 2, response.variable = 3,
                                             n.ahead = n.ahead, atT = atT),
                     name = "") ; p6

p7 <- impres2dplotv2(tmp = get_irf_tvp_unitS(Chain = G000_obj,
                                             impulse.variable = 3, response.variable = 1,
                                             n.ahead = n.ahead, atT = atT),
                     name = "GDP growth") ; p7

p8 <- impres2dplotv2(tmp = get_irf_tvp_unitS(Chain = G000_obj,
                                             impulse.variable = 3, response.variable = 2,
                                             n.ahead = n.ahead, atT = atT),
                     name = "") ; p8
p9 <- impres2dplotv2(tmp = get_irf_tvp_unitS(Chain = G000_obj,
                                             impulse.variable = 3, response.variable = 3,
                                             n.ahead = n.ahead, atT = atT),
                     name = "") ; p9

png(filename='/home/hoanguc3m/Dropbox/WP11/Code/htvpAU/img/IR/IRF_G000_unit_2025Q1.png', width = 1200*0.6, height = 1200*0.6)
#plot_grid(p12, bottom_row, labels = c('', ''), ncol = 1)
par(mar=c(2,1,3,5))
grid.arrange( p1, p2, p3,
              p4, p5, p6,
              p7, p8, p9,
              nrow = 3, ncol = 3, as.table = FALSE)
dev.off()
