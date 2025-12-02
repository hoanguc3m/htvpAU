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

atT <- 110
D_range[atT]
p1 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = T000_obj,
                                            impulse.variable = 1, response.variable = 1,
                                            n.ahead = n.ahead, atT = atT),
                     name = "EPU", cid = "EPU") ; p1

p2 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = T000_obj,
                                            impulse.variable = 1, response.variable = 2,
                                            n.ahead = n.ahead, atT = atT),
                     name = "", cid = "Unemployment") ; p2
p3 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = T000_obj,
                                            impulse.variable = 1, response.variable = 3,
                                            n.ahead = n.ahead, atT = atT),
                     name = "", cid = "GDP growth") ; p3

p4 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = T000_obj,
                                            impulse.variable = 2, response.variable = 1,
                                            n.ahead = n.ahead, atT = atT),
                     name = "Unemployment") ; p4


p5 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = T000_obj,
                                            impulse.variable = 2, response.variable = 2,
                                            n.ahead = n.ahead, atT = atT),
                     name = "") ; p5

p6 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = T000_obj,
                                            impulse.variable = 2, response.variable = 3,
                                            n.ahead = n.ahead, atT = atT),
                     name = "") ; p6

p7 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = T000_obj,
                                            impulse.variable = 3, response.variable = 1,
                                            n.ahead = n.ahead, atT = atT),
                     name = "GDP growth") ; p7

p8 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = T000_obj,
                                            impulse.variable = 3, response.variable = 2,
                                            n.ahead = n.ahead, atT = atT),
                     name = "") ; p8
p9 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = T000_obj,
                                            impulse.variable = 3, response.variable = 3,
                                            n.ahead = n.ahead, atT = atT),
                     name = "") ; p9

png(filename='/home/hoanguc3m/Dropbox/WP11/img/IR/IRF_T000_std_2025Q1.png', width = 1200, height = 1200)
#plot_grid(p12, bottom_row, labels = c('', ''), ncol = 1)
par(mar=c(2,1,3,5))
grid.arrange( p1, p2, p3,
              p4, p5, p6,
              p7, p8, p9,
              nrow = 3, ncol = 3, as.table = FALSE)
dev.off()


##############################################
load("/home/hoanguc3m/Downloads/htvpAU/G000.RData")

#Model G000
atT <- 110
D_range[atT]
p1 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = G000_obj,
                                            impulse.variable = 1, response.variable = 1,
                                            n.ahead = n.ahead, atT = atT),
                     name = "EPU", cid = "EPU") ; p1

p2 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = G000_obj,
                                            impulse.variable = 1, response.variable = 2,
                                            n.ahead = n.ahead, atT = atT),
                     name = "", cid = "Unemployment") ; p2
p3 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = G000_obj,
                                            impulse.variable = 1, response.variable = 3,
                                            n.ahead = n.ahead, atT = atT),
                     name = "", cid = "GDP growth") ; p3

p4 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = G000_obj,
                                            impulse.variable = 2, response.variable = 1,
                                            n.ahead = n.ahead, atT = atT),
                     name = "Unemployment") ; p4


p5 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = G000_obj,
                                            impulse.variable = 2, response.variable = 2,
                                            n.ahead = n.ahead, atT = atT),
                     name = "") ; p5

p6 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = G000_obj,
                                            impulse.variable = 2, response.variable = 3,
                                            n.ahead = n.ahead, atT = atT),
                     name = "") ; p6

p7 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = G000_obj,
                                            impulse.variable = 3, response.variable = 1,
                                            n.ahead = n.ahead, atT = atT),
                     name = "GDP growth") ; p7

p8 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = G000_obj,
                                            impulse.variable = 3, response.variable = 2,
                                            n.ahead = n.ahead, atT = atT),
                     name = "") ; p8
p9 <- impres2dplotv2(tmp = get_irf_tvp_stdS(Chain = G000_obj,
                                            impulse.variable = 3, response.variable = 3,
                                            n.ahead = n.ahead, atT = atT),
                     name = "") ; p9

png(filename='/home/hoanguc3m/Dropbox/WP11/Code/htvpAU/img/IR/IRF_G000_std_2025Q1.png', width = 1200, height = 1200)
#plot_grid(p12, bottom_row, labels = c('', ''), ncol = 1)
par(mar=c(2,1,3,5))
grid.arrange( p1, p2, p3,
              p4, p5, p6,
              p7, p8, p9,
              nrow = 3, ncol = 3, as.table = FALSE)
dev.off()


##############################################
