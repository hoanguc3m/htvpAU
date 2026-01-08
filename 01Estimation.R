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
priors <- get_prior_minnesota(y = y, p = p, intercept=TRUE)

####################
#Tight prior
priors$b0[2] = 0.9
priors$b0[10] = 0.9
priors$b0[18] = 0 # gdp growth

B0 <- t(matrix(priors$b0, ncol = K))
B0[,1] <- c(0.1,0.5,0.0) #/ c(0.2, 0.4, 1)
priors$b0 <- as.numeric(t(B0))

lambda1 <- 0.2
lambda2 <- 0.5
lambda2a <- 0.000001
lambda3 <- 1
lambda4 <- 10

sigmasq <- rep(0, K)
for (ii in 1:K) {
  EstMdl <- arima(y[, ii], order = c(p, 0, 0))
  sigmasq[ii] <- EstMdl$sigma2
}

M <- K * p + 1
Vi <- rep(0, K * M)

for (ii in 1:K) {
    Vi[ii] <- lambda4 * sigmasq[ii]
  }
for (jj in 1:p) { # lag
  for (kk in 1:K) { # B matrix
    for (ii in 1:K) { # B matrix
      indx <- (kk - 1) * K + (jj - 1) * K * K + ii + 1 * K
      if (ii == kk) {
        Vi[indx] <- lambda1/jj^2 # Own lag
      } else {
        if (ii > 1) {
          Vi[indx] <- (lambda1 * lambda2)/(jj^2) * (sigmasq[ii]/sigmasq[kk]) # Cross lag of EPU to U or GDP
        }
        if (ii == 1) {
          Vi[indx] <- (lambda1 * lambda2a)/(jj^2) * (sigmasq[ii]/sigmasq[kk]) #Cross lag of U or GDP to EPU
        }
      }
    }
  }
}
V_b_prior = as.numeric(t(matrix(Vi, nrow = K)))

priors$V_b_prior <- V_b_prior
priors$V_a0_prior <- rep(0.01,length(priors$a0))
priors$hyper_ab <- list(hyper_a = rep(0.0001, 3),
                        hyper_b = c(c(1e-4,1e-4,1e-10,1e-10,1e-4,1e-10,1e-10),
                                    rep(0.0001,21)) )

#################################
inits <- list(samples = 100000, burnin = 20000, thin = 20)
RhpcBLASctl::blas_set_num_threads(2)
setwd("/home/hoanguc3m/Downloads/htvpAU")
{
  inits$is_tv = c(0,0,0); G000_obj <- fitTVPGaussSV_ASIS(y, y0, p, priors, inits)
  save(G000_obj, file = "G000.RData")

  # inits$is_tv = c(0,0,0); G000_NCP <- fitTVPGaussSV_NCP(y, y0, p, priors, inits)
  # save(G000_NCP, file = "G000_NCP.RData")

  inits$is_tv = c(0,0,1); G001_obj <- fitTVPGaussSV_ASIS(y, y0, p, priors, inits)
  save(G001_obj, file = "G001.RData")


  inits$is_tv = c(0,1,0); G010_obj <- fitTVPGaussSV_ASIS(y, y0, p, priors, inits)
  save(G010_obj, file = "G010.RData")

  inits$is_tv = c(1,0,0); G100_obj <- fitTVPGaussSV_ASIS(y, y0, p, priors, inits)
  save(G100_obj, file = "G100.RData")

  inits$is_tv = c(1,1,0); G110_obj <- fitTVPGaussSV_ASIS(y, y0, p, priors, inits)
  save(G110_obj, file = "G110.RData")

  inits$is_tv = c(0,1,1); G011_obj <- fitTVPGaussSV_ASIS(y, y0, p, priors, inits)
  save(G011_obj, file = "G011.RData")

  inits$is_tv = c(1,0,1); G101_obj <- fitTVPGaussSV_ASIS(y, y0, p, priors, inits)
  save(G101_obj, file = "G101.RData")
  #inits$samples <- 20000
  # inits$samples <- 20000
  # inits$thin <- 1
  inits$is_tv = c(1,1,1); G111_obj <- fitTVPGaussSV_ASIS(y, y0, p, priors, inits)
  save(G111_obj, file = "G111.RData")

  # inits$is_tv = c(1,1,1); G111_NCP <- fitTVPGaussSV_NCP(y, y0, p, priors, inits)
  # save(G111_NCP, file = "G111_NCP.RData")

}
####################################################################
{
  inits$is_tv = c(0,0,0); T000_obj <- fitTVPStudentSV_ASIS(y, y0, p, priors, inits)
  save(T000_obj, file = "T000.RData")

  inits$is_tv = c(0,0,1); T001_obj <- fitTVPStudentSV_ASIS(y, y0, p, priors, inits)
  save(T001_obj, file = "T001.RData")

  inits$is_tv = c(0,1,0); T010_obj <- fitTVPStudentSV_ASIS(y, y0, p, priors, inits)
  save(T010_obj, file = "T010.RData")

  inits$is_tv = c(1,0,0); T100_obj <- fitTVPStudentSV_ASIS(y, y0, p, priors, inits)
  save(T100_obj, file = "T100.RData")

  inits$is_tv = c(1,1,0); T010_obj <- fitTVPStudentSV_ASIS(y, y0, p, priors, inits)
  save(T010_obj, file = "T010.RData")

  inits$is_tv = c(0,1,1); T011_obj <- fitTVPStudentSV_ASIS(y, y0, p, priors, inits)
  save(T011_obj, file = "T011.RData")

  inits$is_tv = c(1,0,1); T101_obj <- fitTVPStudentSV_ASIS(y, y0, p, priors, inits)
  save(T101_obj, file = "T101.RData")

  inits$is_tv = c(1,1,1); T111_obj <- fitTVPStudentSV_ASIS(y, y0, p, priors, inits)
  save(T111_obj, file = "T111.RData")
}

###############################################
setwd("/home/hoanguc3m/Downloads/htvpAU")
load("T111.RData")
load("T000.RData")

####################################

library(ggplot2)
library(gridExtra)
library(cowplot)
library(ggthemr)
ggthemr('light')

Time <- seq(as.Date("1997/06/01"), as.Date("2025/09/01"), "quarter")
recessions.df = read.table(textConnection(
  "Peak, Trough
  2001-03-01, 2001-11-01
  2007-12-01, 2009-06-01
  2020-03-01, 2020-05-01"), sep=',',
  colClasses=c('Date', 'Date'), header=TRUE)
data_nu <- data.frame(Time = Time, EPU = dataraw$GLOBAL_EPU,
                      U = dataraw$U,
                      GDP = dataraw$GROWTH)

p1 <- ggplot(data_nu, aes(x = Time, y = EPU)) +
  geom_line(color = "#e3625b") + xlab("EPU") + ylab("") + ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) +
  geom_rect(data=recessions.df, inherit.aes=F, aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill='darkgray', alpha=0.5) + theme_bw()

p2 <- ggplot(data_nu, aes(x = Time, y = U)) +
  geom_line(color = "#e3625b") + xlab("Unemployment") + ylab("") +
  geom_rect(data=recessions.df, inherit.aes=F, aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill='darkgray', alpha=0.5) + theme_bw()

p3 <- ggplot(data_nu, aes(x = Time, y = GDP)) +
  geom_line(color = "#e3625b") + xlab("GDP growth") + ylab("") +
  geom_rect(data=recessions.df, inherit.aes=F, aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill='darkgray', alpha=0.5) + theme_bw()

setwd("/home/hoanguc3m/Dropbox/WP11/Code/htvpAU/")
pdf(file='img/data.pdf', width = 8.1, height = 5.4)
plot_grid(p1, p2, p3, ncol = 1, align = "v")
#grid.arrange(p1, p2, p3, nrow = 3, ncol = 1)
dev.off()



library(gridExtra)
library(ggthemr)
# ggthemr('flat')
ggthemr('solarized')


Nu_mat111 <- T111_obj$store_nu
Nu_mat000 <- T000_obj$store_nu
ndraws <- nrow(Nu_mat000)
varname <- c("EPU", "Unemployemnt", "GDP growth")
gg_nu_mat <- function(i){
  data_nu <- data.frame(nu = c(Nu_mat000[,i], Nu_mat111[,i]), Models = c(rep("T000", ndraws), rep("T111", ndraws)))
  # max_nu <- max(c(Nu_mat000[,i], Nu_mat111[,i]))
  max_nu <- 75
  data_ann <- data.frame( x = 0.7*min(data_nu[,1]) + 0.3*max(data_nu[,1]),
                          y = c(max(density(data_nu[,1])$y) *0.5 + 0.2 * sd(density(data_nu[,1])$y),
                                max(density(data_nu[,1])$y) *0.5 - 0.2 * sd(density(data_nu[,1])$y)),
                          label = c(paste("T000 : Mean =", round(mean(c(Nu_mat000[,i])),1), ", Median = ", round(median(c(Nu_mat000[,i])),1)),
                                    paste("T111 : Mean =", round(mean(c(Nu_mat111[,i])),1), ", Median = ", round(median(c(Nu_mat111[,i])),1)) ),
                          Models = c("T000","T111"))

  ggplot(data_nu) +
    geom_density(aes(x=nu, colour = Models, fill = Models, linetype = Models), alpha = 0.25, adjust = 2, size=1.1 ) +
    scale_colour_ggthemr_d() + ylab(expression(nu)) + xlab(varname[i])  +
    xlim(4, max_nu) + theme_light() + scale_linetype_manual(values=c("dotted", "longdash")) + theme(legend.position="bottom") +
    #round(mean(c(Nu_mat000[,i])),1)
    # annotate(geom="text", x=50, y=0.05, label=paste("M =", round(mean(c(Nu_mat111[,i])),1), " Mdn = ", round(median(c(Nu_mat111[,i])),1))  , color="red") +
    # annotate(geom="text", x=50, y=0.04, label=paste("M =", round(mean(c(Nu_mat000[,i])),1), " Mdn = ", round(median(c(Nu_mat000[,i])),1))  , color="blue") +
    geom_text(data=data_ann, aes( x=x, y=y, color=Models, label=label), hjust = 0,
              fontface="bold", size = 2.5 )



}

p1 <- gg_nu_mat(1)
p2 <- gg_nu_mat(2)
p3 <- gg_nu_mat(3)


pdf(file='img/postNu.pdf', width = 12, height = 4)
par(mar=c(2,5,3,1))
grid.arrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()
####################################################################
Sigma_h111 <- T111_obj$store_Sigh
Sigma_h000 <- T000_obj$store_Sigh
gg_sigma_mat <- function(i){
  data_nu <- data.frame(nu = c(Sigma_h000[,i], Sigma_h111[,i]), Models = c(rep("T000", ndraws), rep("T111", ndraws)))
  max_nu <- max(c(Sigma_h000[,i], Sigma_h111[,i]))
  data_ann <- data.frame( x = 0.5*min(data_nu[,1]) + 0.5*max(data_nu[,1]),
                          y = c(max(density(data_nu[,1])$y) *0.5 + 0.2 * sd(density(data_nu[,1])$y),
                                max(density(data_nu[,1])$y) *0.5 - 0.2 * sd(density(data_nu[,1])$y)),
                          label = c(paste("T000 : Mean =", round(mean(c(Sigma_h000[,i])),2), ", Median = ", round(median(c(Sigma_h000[,i])),2)),
                                    paste("T111 : Mean =", round(mean(c(Sigma_h111[,i])),2), ", Median = ", round(median(c(Sigma_h111[,i])),2)) ),
                          Models = c("T000","T111"))

  ggplot(data_nu) +
    geom_density(aes(x=nu, colour = Models, fill = Models, linetype = Models), alpha = 0.25, adjust = 1.5, size=1.1 ) +
    scale_colour_ggthemr_d() +
    ylab(expression(sigma[h]^2)) + xlab(varname[i])  +
    theme_light() + scale_linetype_manual(values=c("dotted", "longdash")) + theme(legend.position="bottom") +
    geom_text(data=data_ann, aes( x=x, y=y, color=Models, label=label), hjust = 0,
              fontface="bold", size = 2.5 )

}
p1 <- gg_sigma_mat(1)
p2 <- gg_sigma_mat(2)
p3 <- gg_sigma_mat(3)

pdf(file='img/postSigmah.pdf', width = 12, height = 4)
par(mar=c(2,5,3,1))
grid.arrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()

####################################################################
Sigma_ab111 <- cbind(T111_obj$store_Sigbeta, T111_obj$store_Sigalp)
varname <- c("b[1]", "B1[1,1]", "B1[1,2]", "B1[1,3]", "B2[1,1]", "B2[1,2]", "B2[1,3]",
             "b[2]", "B1[2,1]", "B1[2,2]", "B1[2,3]", "B2[2,1]", "B2[2,2]", "B2[2,3]",
             "b[3]", "B1[3,1]", "B1[3,2]", "B1[3,3]", "B2[3,1]", "B2[3,2]", "B2[3,3]",
             "A[2,1]", "A[3,1]", "A[3,2]")

gg_sigmaAB_mat <- function(i){
  data_nu <- data.frame(nu = c((Sigma_ab111[,i])) , Models = c(rep("T111", ndraws)))
  max_nu <- (max(c(Sigma_ab111[,i])))
  ggplot(data_nu) +
    geom_density(aes(x=nu, colour = Models, fill = Models), alpha = 0.5) +
    scale_colour_ggthemr_d() +
    xlim(0,max_nu) + #geom_histogram(position = "identity", alpha = 0.8, bins = 30) +
    xlab(bquote(.(varname[i]))) + ylab(expression(sigma^2)) + theme_bw() + theme(legend.position="bottom")
}
l <- list()
for (i in c(1:24)) l[[i]] <- gg_sigmaAB_mat(i)


pdf(file='img/postSigmaAB.pdf', width = 12, height = 3*8)
par(mar=c(2,5,3,1))
plot_grid(l[[1]], l[[8]], l[[15]],
          l[[2]], l[[9]], l[[16]],
          l[[3]], l[[10]], l[[17]],
          l[[4]], l[[11]], l[[18]],
          l[[5]], l[[12]], l[[19]],
          l[[6]], l[[13]], l[[20]],
          l[[7]], l[[14]], l[[21]],
          l[[22]], l[[23]], l[[24]],
          ncol = 3, byrow = TRUE)
dev.off()

####################################################################
ab111_mean <- cbind( apply(get_post(T111_obj,element = "beta"), MARGIN = c(2,3), FUN = mean),
                     apply(get_post(T111_obj,element = "alpha"), MARGIN = c(2,3), FUN = mean))
ab111_q10 <- cbind( apply(get_post(T111_obj,element = "beta"), MARGIN = c(2,3), FUN = quantile, probs = c(0.1)),
                    apply(get_post(T111_obj,element = "alpha"), MARGIN = c(2,3), FUN = quantile, probs = c(0.1)))
ab111_q90 <- cbind( apply(get_post(T111_obj,element = "beta"), MARGIN = c(2,3), FUN = quantile, probs = c(0.9)),
                    apply(get_post(T111_obj,element = "alpha"), MARGIN = c(2,3), FUN = quantile, probs = c(0.9)))

gg_TVPAB_mat <- function(i){
  Time <- tail(seq(as.Date("1997/06/01"), as.Date("2025/09/01"), "quarter"), nrow(ab111_mean))

  data_nu <- data.frame(Time = Time, AB_mean = ab111_mean[,i])

  data_nu1 <- data.frame(Time = c(Time, rev(Time)) ,
                         AB_UL = c(ab111_q10[,i], rev(ab111_q90[,i])))
  x_lab_name <- bquote(.(varname[i]))
  ggplot() + geom_line(data=data_nu, mapping=aes(x=Time, y=AB_mean), color = "#ff0c00") + scale_colour_ggthemr_d() +
    geom_polygon(data=data_nu1, mapping=aes(x=Time, y=AB_UL), fill = "#5ba6d6", alpha = 0.5) +
    xlab(varname[i]) + ylab("") + theme_bw() + theme(legend.position="bottom") + theme(plot.title = element_text(hjust = 0.5)) +
    geom_rect(data=recessions.df, inherit.aes=F, aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill='darkgray', alpha=0.5)

}
varname <- parse(text = c("Intercept", "e[t-1]", "u[t-1]", "g[t-1]", "e[t-2]", "u[t-2]", "g[t-2]",
                          "Intercept", "e[t-1]", "u[t-1]", "g[t-1]", "e[t-2]", "u[t-2]", "g[t-2]",
                          "Intercept", "e[t-1]", "u[t-1]", "g[t-1]", "e[t-2]", "u[t-2]", "g[t-2]",
                          "a[21,t]", "a[31,t]", "a[32,t]"))
l <- list()
for (i in c(1:24)) l[[i]] <- gg_TVPAB_mat(i)

pdf(file='img/postAB_eq1.pdf', width = 12*0.7, height = 3*4*0.7)
par(mar=c(2,5,3,1))
plot_grid(l[[2]], l[[5]],
          l[[3]], l[[6]],
          l[[4]], l[[7]],
          l[[1]], NULL,
          ncol = 2, byrow = TRUE)
dev.off()

pdf(file='img/postAB_eq2.pdf', width = 12*0.7, height = 3*4*0.7)
par(mar=c(2,5,3,1))
plot_grid(l[[9]], l[[12]],
          l[[10]], l[[13]],
          l[[11]], l[[14]],
          l[[8]], l[[22]],
          ncol = 2, byrow = TRUE)
dev.off()

pdf(file='img/postAB_eq3.pdf', width = 12*0.7, height = 3*5*0.7)
par(mar=c(2,5,3,1))
plot_grid(l[[16]], l[[19]],
          l[[17]], l[[20]],
          l[[18]], l[[21]],
          l[[15]], NULL,
          l[[23]], l[[24]],
          ncol = 2, byrow = TRUE)
dev.off()

pdf(file='img/postAB.pdf', width = 12, height = 3*8)
par(mar=c(2,5,3,1))
plot_grid(l[[1]], l[[8]], l[[15]],
          l[[2]], l[[9]], l[[16]],
          l[[3]], l[[10]], l[[17]],
          l[[4]], l[[11]], l[[18]],
          l[[5]], l[[12]], l[[19]],
          l[[6]], l[[13]], l[[20]],
          l[[7]], l[[14]], l[[21]],
          l[[22]], l[[23]], l[[24]],
          ncol = 3, byrow = TRUE)
dev.off()

pdf(file='img/postAB_long.pdf', height = 12, width = 3*7)
par(mar=c(2,5,3,1))
plot_grid(l[[1]], l[[8]], l[[15]],
          l[[2]], l[[9]], l[[16]],
          l[[3]], l[[10]], l[[17]],
          l[[4]], l[[11]], l[[18]],
          l[[5]], l[[12]], l[[19]],
          l[[6]], l[[13]], l[[20]],
          l[[7]], l[[14]], l[[21]],
          l[[22]], l[[23]], l[[24]],
          nrow = 3, byrow = FALSE)
dev.off()
#####################################################################
{
  beta_reduce_post <- get_post(T111_obj,element = "betareduce")
  rho_reduce_post <- get_post(T111_obj,element = "rho")
}
{

  ab111_mean <- cbind( apply(beta_reduce_post, MARGIN = c(2,3), FUN = mean) )
  ab111_q10 <- cbind( apply(beta_reduce_post, MARGIN = c(2,3), FUN = quantile, probs = c(0.1)) )
  ab111_q90 <- cbind( apply(beta_reduce_post, MARGIN = c(2,3), FUN = quantile, probs = c(0.9)) )

  gg_TVP_redu_AB_mat <- function(i){
    Time <- tail(seq(as.Date("1997/06/01"), as.Date("2025/09/01"), "quarter"), nrow(ab111_mean))

    data_nu <- data.frame(Time = Time, AB_mean = ab111_mean[,i])

    data_nu1 <- data.frame(Time = c(Time, rev(Time)) ,
                           AB_UL = c(ab111_q10[,i], rev(ab111_q90[,i])))
    x_lab_name <- bquote(.(varname[i]))
    ggplot() + geom_line(data=data_nu, mapping=aes(x=Time, y=AB_mean), color = "#ff0c00") + scale_colour_ggthemr_d() +
      geom_polygon(data=data_nu1, mapping=aes(x=Time, y=AB_UL), fill = "#5ba6d6", alpha = 0.5) +
      xlab(varname[i]) + ylab("") + theme_bw() + theme(legend.position="bottom") + theme(plot.title = element_text(hjust = 0.5)) +
      geom_rect(data=recessions.df, inherit.aes=F, aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill='darkgray', alpha=0.5)

  }
  varname <- parse(text = c("Intercept", "e[t-1]", "u[t-1]", "g[t-1]", "e[t-2]", "u[t-2]", "g[t-2]",
                            "Intercept", "e[t-1]", "u[t-1]", "g[t-1]", "e[t-2]", "u[t-2]", "g[t-2]",
                            "Intercept", "e[t-1]", "u[t-1]", "g[t-1]", "e[t-2]", "u[t-2]", "g[t-2]"))
  l <- list()
  for (i in c(1:21)) l[[i]] <- gg_TVP_redu_AB_mat(i)



  pdf(file='img/post_red_AB_eq1.pdf', width = 12*0.7, height = 3*4*0.7)
  par(mar=c(2,5,3,1))
  plot_grid(l[[2]], l[[5]],
            l[[3]], l[[6]],
            l[[4]], l[[7]],
            l[[1]], NULL,
            ncol = 2, byrow = TRUE)
  dev.off()

  pdf(file='img/post_red_AB_eq2.pdf', width = 12*0.7, height = 3*4*0.7)
  par(mar=c(2,5,3,1))
  plot_grid(l[[9]], l[[12]],
            l[[10]], l[[13]],
            l[[11]], l[[14]],
            l[[8]], NULL,
            ncol = 2, byrow = TRUE)
  dev.off()

  pdf(file='img/post_red_AB_eq3.pdf', width = 12*0.7, height = 3*4*0.7)
  par(mar=c(2,5,3,1))
  plot_grid(l[[16]], l[[19]],
            l[[17]], l[[20]],
            l[[18]], l[[21]],
            l[[15]],  NULL,
            ncol = 2, byrow = TRUE)
  dev.off()


  rho111_mean <- ( apply(rho_reduce_post, MARGIN = c(2,3), FUN = mean))
  rho111_q10 <- ( apply(rho_reduce_post, MARGIN = c(2,3), FUN = quantile, probs = c(0.1)))
  rho111_q90 <- ( apply(rho_reduce_post, MARGIN = c(2,3), FUN = quantile, probs = c(0.9)) )
  varname <- parse(text = c("rho[paste(2, ",", 1)]", "rho[paste(3, ",", 1)]", "rho[paste(3, ",", 2)]"))
  gg_H_mat <- function(i){
    Time <- tail(seq(as.Date("1997/06/01"), as.Date("2025/09/01"), "quarter"), nrow(rho111_mean))

    data_nu <- data.frame(Time = Time, H_mean = rho111_mean[,i])

    data_nu1 <- data.frame(Time = c(Time, rev(Time)) ,
                           H_UL = c(rho111_q10[,i], rev(rho111_q90[,i])))
    miny <- min(rho111_q10[,i]) - 0.1
    maxy <- max(rho111_q90[,i]) + 0.1

    ggplot() + geom_line(data=data_nu, mapping=aes(x=Time, y=H_mean), color = "#ff0c00") + ylim(c(miny, maxy)) +
      geom_polygon(data=data_nu1, mapping=aes(x=Time, y=H_UL), fill = "#5ba6d6", alpha = 0.5) +
      xlab(varname[i]) + ylab("") + theme_bw() + theme(legend.position="bottom") + theme(plot.title = element_text(hjust = 0.5)) +
      geom_rect(data=recessions.df, inherit.aes=F, aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill='darkgray', alpha=0.5)

  }

  l <- list()
  for (i in c(1:3)) l[[i]] <- gg_H_mat(i)




  pdf(file='img/post_redu_rho.pdf', width = 12*0.7, height = 3*0.7)
  par(mar=c(2,5,3,1))
  plot_grid(l[[1]], l[[2]], l[[3]],
            ncol = 3, byrow = TRUE)
  dev.off()
}


####################################################################
load("T010.RData")

h000_mean <- ( apply((get_post(T010_obj,element = "h")), MARGIN = c(2,3), FUN = mean))
h000_q10 <- ( apply((get_post(T010_obj,element = "h")), MARGIN = c(2,3), FUN = quantile, probs = c(0.1)))
h000_q90 <- ( apply((get_post(T010_obj,element = "h")), MARGIN = c(2,3), FUN = quantile, probs = c(0.9)) )

# h000_mean <- ( apply(exp(get_post(T010_obj,element = "h")), MARGIN = c(2,3), FUN = mean))
# h000_q10 <- ( apply(exp(get_post(T010_obj,element = "h")), MARGIN = c(2,3), FUN = quantile, probs = c(0.1)))
# h000_q90 <- ( apply(exp(get_post(T010_obj,element = "h")), MARGIN = c(2,3), FUN = quantile, probs = c(0.9)) )

varname <- c("EPU", "Unemployment", "GDP growth")
gg_H_mat <- function(i){
  Time <- tail(seq(as.Date("1997/06/01"), as.Date("2025/09/01"), "quarter"), nrow(h000_mean))

  data_nu <- data.frame(Time = Time, H_mean = h000_mean[,i])

  data_nu1 <- data.frame(Time = c(Time, rev(Time)) ,
                         H_UL = c(h000_q10[,i], rev(h000_q90[,i])))
  miny <- min(h000_q10[,i]) #- 0.3
  maxy <- max(h000_q90[,i]) + 0.3

  ggplot() + geom_line(data=data_nu, mapping=aes(x=Time, y=H_mean), color = "#ff0c00") + ylim(c(miny, maxy)) +
    geom_polygon(data=data_nu1, mapping=aes(x=Time, y=H_UL), fill = "#5ba6d6", alpha = 0.5) +
    xlab(bquote(.(varname[i]))) + ylab("") + theme_bw() + theme(legend.position="bottom") + theme(plot.title = element_text(hjust = 0.5)) +
    geom_rect(data=recessions.df, inherit.aes=F, aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill='darkgray', alpha=0.5)

}

l <- list()
for (i in c(1:3)) l[[i]] <- gg_H_mat(i)


pdf(file='img/postHT010.pdf', width = 9, height = 6)
plot_grid(l[[1]], l[[2]], l[[3]],
          ncol = 1, align = "v")
#grid.arrange(p1, p2, p3, nrow = 3, ncol = 1)
dev.off()

####################################################################
load("T010.RData")

h000_mean <- ( apply(exp(get_post(T010_obj,element = "h")), MARGIN = c(2,3), FUN = mean))
h000_q10 <- ( apply(exp(get_post(T010_obj,element = "h")), MARGIN = c(2,3), FUN = quantile, probs = c(0.1)))
h000_q90 <- ( apply(exp(get_post(T010_obj,element = "h")), MARGIN = c(2,3), FUN = quantile, probs = c(0.9)) )

# h000_mean <- ( apply(exp(get_post(T010_obj,element = "h")), MARGIN = c(2,3), FUN = mean))
# h000_q10 <- ( apply(exp(get_post(T010_obj,element = "h")), MARGIN = c(2,3), FUN = quantile, probs = c(0.1)))
# h000_q90 <- ( apply(exp(get_post(T010_obj,element = "h")), MARGIN = c(2,3), FUN = quantile, probs = c(0.9)) )

varname <- c("EPU", "Unemployment", "GDP growth")
gg_H_mat <- function(i){
  Time <- tail(seq(as.Date("1997/06/01"), as.Date("2025/09/01"), "quarter"), nrow(h000_mean))

  data_nu <- data.frame(Time = Time, H_mean = h000_mean[,i])

  data_nu1 <- data.frame(Time = c(Time, rev(Time)) ,
                         H_UL = c(h000_q10[,i], rev(h000_q90[,i])))
  miny <- min(h000_q10[,i]) #- 0.3
  maxy <- max(h000_q90[,i]) + 0.3

  ggplot() + geom_line(data=data_nu, mapping=aes(x=Time, y=H_mean), color = "#ff0c00") + ylim(c(miny, maxy)) +
    geom_polygon(data=data_nu1, mapping=aes(x=Time, y=H_UL), fill = "#5ba6d6", alpha = 0.5) +
    xlab(bquote(.(varname[i]))) + ylab("") + theme_bw() + theme(legend.position="bottom") + theme(plot.title = element_text(hjust = 0.5)) +
    geom_rect(data=recessions.df, inherit.aes=F, aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill='darkgray', alpha=0.5)

}

l <- list()
for (i in c(1:3)) l[[i]] <- gg_H_mat(i)

library(cowplot)

pdf(file='img/poststructureHT010.pdf', width = 9, height = 6)
plot_grid(l[[1]], l[[2]], l[[3]],
          ncol = 1, align = "v")
#grid.arrange(p1, p2, p3, nrow = 3, ncol = 1)
dev.off()

####################################################################
load("T010.RData")

h000_mean <- ( apply(exp(get_post(T010_obj,element = "h")), MARGIN = c(2,3), FUN = mean))
h000_q10 <- ( apply(exp(get_post(T010_obj,element = "h")), MARGIN = c(2,3), FUN = quantile, probs = c(0.1)))
h000_q90 <- ( apply(exp(get_post(T010_obj,element = "h")), MARGIN = c(2,3), FUN = quantile, probs = c(0.9)) )

# h000_mean <- ( apply(exp(get_post(T010_obj,element = "h")), MARGIN = c(2,3), FUN = mean))
# h000_q10 <- ( apply(exp(get_post(T010_obj,element = "h")), MARGIN = c(2,3), FUN = quantile, probs = c(0.1)))
# h000_q90 <- ( apply(exp(get_post(T010_obj,element = "h")), MARGIN = c(2,3), FUN = quantile, probs = c(0.9)) )

varname <- c("EPU", "Unemployment", "GDP growth")
gg_H_mat <- function(i){
  Time <- tail(seq(as.Date("1997/06/01"), as.Date("2025/09/01"), "quarter"), nrow(h000_mean))

  data_nu <- data.frame(Time = Time, H_mean = h000_mean[,i])

  data_nu1 <- data.frame(Time = c(Time, rev(Time)) ,
                         H_UL = c(h000_q10[,i], rev(h000_q90[,i])))
  miny <- min(h000_q10[,i]) #- 0.3
  maxy <- max(h000_q90[,i]) + 0.3

  ggplot() + geom_line(data=data_nu, mapping=aes(x=Time, y=H_mean), color = "#ff0c00") + ylim(c(miny, maxy)) +
    geom_polygon(data=data_nu1, mapping=aes(x=Time, y=H_UL), fill = "#5ba6d6", alpha = 0.5) +
    xlab(bquote(.(varname[i]))) + ylab("") + theme_bw() + theme(legend.position="bottom") + theme(plot.title = element_text(hjust = 0.5)) +
    geom_rect(data=recessions.df, inherit.aes=F, aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill='darkgray', alpha=0.5)

}

l <- list()
for (i in c(1:3)) l[[i]] <- gg_H_mat(i)


pdf(file='img/postexpHT001.pdf', width = 9, height = 6)
plot_grid(l[[1]], l[[2]], l[[3]],
          ncol = 1, align = "v")
#grid.arrange(p1, p2, p3, nrow = 3, ncol = 1)
dev.off()
###########################################################


get_ht <- function(Chain, atT = NULL, n.ahead = 0, structure = TRUE){
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

  H.chol <- array(diag(K), dim = c(K,K, n.ahead+1))
  beta <- get_post(Chain,element = "beta")
  alp <- get_post(Chain,element = "alpha")

  beta <- beta[,,-c( seq(from = 1, to = k_beta, by = k_beta_div_K) )] # Remove constant


  out_all <- matrix(NA, nrow = K, ncol = ndraws)

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

    s = 1
    sigtemp <- sig[samp,c(1:K)]

    if (structure){
      H.chol[,,s] <- diag(sigtemp)
    } else {
      H.chol[,,s] <- A_inv %*% diag(sigtemp)
    }

    out_all[,j] <- diag(H.chol[,,s])

  }
  return(list(H.chol = out_all )  )
}

T010_diagH_stdS <- parallel::mclapply(c(1:nrow(T010_obj$data$y)),
                                      FUN = function(atT) {
                                        get_ht(Chain = T010_obj,
                                               atT = atT, structure = TRUE)$H.chol
                                      },
                                      mc.cores = 16)
T010_diagH_stdS_structure <- array(NA, dim = dim(get_post(T010_obj,element = "h")))
atT <- nrow(y)
for (t in c(1: (atT) )){
  T010_diagH_stdS_structure[,t,] <- t(T010_diagH_stdS[[t]])
}

h000_mean <- ( apply(T010_diagH_stdS_structure, MARGIN = c(2,3), FUN = mean, na.rm = T))
h000_q10 <- ( apply(T010_diagH_stdS_structure, MARGIN = c(2,3), FUN = quantile, probs = c(0.1), na.rm = T))
h000_q90 <- ( apply(T010_diagH_stdS_structure, MARGIN = c(2,3), FUN = quantile, probs = c(0.9), na.rm = T) )

varname <- c("EPU", "Unemployment", "GDP growth")
gg_H_mat <- function(i){
  Time <- tail(seq(as.Date("1997/06/01"), as.Date("2025/09/01"), "quarter"), nrow(h000_mean))

  data_nu <- data.frame(Time = Time, H_mean = h000_mean[,i])

  data_nu1 <- data.frame(Time = c(Time, rev(Time)) ,
                         H_UL = c(h000_q10[,i], rev(h000_q90[,i])))
  miny <- 0 # min(h000_q10[,i]) #- 0.3
  maxy <- max(h000_q90[,i]) + 0.2
  if(maxy < 0.5) maxy <- 0.35
  # maxy <- 2.5

  ggplot() + geom_line(data=data_nu, mapping=aes(x=Time, y=H_mean), color = "#ff0c00") + ylim(c(miny, maxy)) +
    geom_polygon(data=data_nu1, mapping=aes(x=Time, y=H_UL), fill = "#5ba6d6", alpha = 0.5) +
    xlab(bquote(.(varname[i]))) + ylab("") + theme_bw() + theme(legend.position="bottom") + theme(plot.title = element_text(hjust = 0.5)) +
    geom_rect(data=recessions.df, inherit.aes=F, aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill='darkgray', alpha=0.5)

}

l <- list()
for (i in c(1:3)) l[[i]] <- gg_H_mat(i)

pdf(file='img/postexpHT001div2.pdf', width = 9, height = 6)
plot_grid(l[[1]], l[[2]], l[[3]],
          ncol = 1, align = "v")
#grid.arrange(p1, p2, p3, nrow = 3, ncol = 1)
dev.off()

