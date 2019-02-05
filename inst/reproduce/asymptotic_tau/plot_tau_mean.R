# First, some initialization
# remove all objects from R environment
rm(list = ls())
set.seed(17)
if (!require(CoupledPIMH)){
  library(devtools)
  devtools::document()
}

library(doRNG)
setmytheme()
set.seed(17)


### this script takes around 30 mins on a 90 core server
s_arr <- seq(0.01,5,length.out=100)

# compare actual and predicted meeting times
# for predicted values we need to do a series of integrations over the range of tau appearing in the histograms
# define the conditional geometric probability
p_z <- function(z){return(pnorm(-z/sqrt(s2)) + exp(-z+s2/2)*pnorm((z-s2)/sqrt(s2)))}

I_mean_tau <- function(z){P_z <- (1/p_z(z)) * dnorm(z,0,sd=sqrt(s2))}
I_var_tau <- function(z){P_z <- ((1-p_z(z))/p_z(z)^2) * dnorm(z,0,sd=sqrt(s2))}

tail_arr <- c(1)
p_mat <- matrix(NA,length(tail_arr),length(s_arr))
tau_mean_arr <- rep(NA,length(s_arr))
tau_var_arr <- rep(NA,length(s_arr))
# ok for each estimate of the variance of the log-likelihood we need to estimate the respective tau probabilities
for(i in 1:length(s_arr)){
  s2 <- s_arr[i]^2
  
  lower = -50*sqrt(s2)
  upper = 50*sqrt(s2)

  for(j in 1:length(tail_arr)){
    tail_val <- tail_arr[j]
    I_z <- function(z){P_z <- (1-(1-p_z(z))^tail_val) * dnorm(z,0,sd=sqrt(s2))}
    p_mat[j,i] <- integrate(I_z, lower, upper)$value
  }
  tau_mean_arr[i] <- integrate(I_mean_tau, lower, upper)$value
  tau_var_arr[i] <- integrate(I_var_tau, lower, upper)$value
}

matplot(s_arr,t(p_mat),pch=1:length(tail_arr),type='l')
plot(s_arr,tau_var_arr,type='l',col='orange',log='y')
plot(s_arr,tau_mean_arr,log='',type='l',col='orange')

# R <- 1e4
# zmeans <- rep(NA,R)
# p_z <- function(z){return(pnorm(-z/sqrt(s2)) + exp(-z+s2/2)*pnorm((z-s2)/sqrt(s2)))}
# I_mean_tau <- function(z){P_z <- (1/p_z(z)) * dnorm(z,0,sd=sqrt(s2))}
# s2 <- 25
# zmeans <- foreach(r = 1:R,.combine = rbind) %dorng% {
#   #print(r)
#   z <- rnorm(1e7,0,sd=sqrt(s2))
#   mean(1/p_z(z))
# }


p_mat_df <- data.frame(s=s_arr,p=t(p_mat))
g1<-ggplot(data = p_mat_df) +
   geom_line(aes(x=s,y=p))+
   #geom_errorbar(data=plot_df, aes(x=tau_val,y=value,fill=N,ymin=l,ymax=u),position='dodge') +
   xlab(TeX('$\\sigma$'))+
   scale_x_continuous(breaks=seq(0,10,1)) +
   #ylab()
  ylab('')+
  ggtitle(TeX('$P\\[\\tau = 1\\]$'))+
  theme(plot.title = element_text(size = 20),
        #plot.margin = margin(0, 0, 0, 0),
        axis.title.y = element_blank())#, face = "bold"))

g1

p_mat_df <- data.frame(s=s_arr,p=tau_mean_arr)
g2<-ggplot(data = p_mat_df) +
  geom_line(aes(x=s,y=p))+
  #geom_errorbar(data=plot_df, aes(x=tau_val,y=value,fill=N,ymin=l,ymax=u),position='dodge') +
  scale_x_continuous(breaks=seq(0,10,1)) +
  xlab(TeX('$\\sigma$'))+
  ylab('')+
  ggtitle(TeX('$E\\[\\tau\\]$'))+
  theme(plot.title = element_text(size = 20),
        axis.title.y = element_blank())
g2



library(egg)      
g <- ggarrange(g1,g2,ncol=2)
g
# g <- cowplot::plot_grid(g1, g2, labels = c("",""))
# g
name <- ('ptau_summary.pdf')
ggsave(name,plot=g,width=7,height=4)

# #name <- ('ptau_summary.pdf')
# #ggsave(name,plot=g,width=7,height=5)
# 
# all_df <- data.frame(s=s_arr,p=t(p_mat),mean=tau_mean_arr)
# all_df_m <- reshape2::melt(all_df,id='s')
# all_df_m$variable <- as.factor(all_df_m$variable)
# levels(all_df_m$variable) <- c('P[tau=1]','E[tau]')
# 
# g <- ggplot(data=all_df_m) + 
#   geom_line(aes(x=s,y=value))+
#   facet_wrap(~variable,scales='free')+
#   xlab(TeX('$\\sigma$'))+
#   ylab('')
# g
# name <- ('ptau_summary.pdf')
# ggsave(name,plot=g,width=7,height=5)
# 
# #plot(s2_arr,tau_var_arr^0.5,log='',type='l',col='orange')
# 

# plot(s2_arr,tau_var_arr,log='y',type='l')
# matplot(s2_arr,tau_mean_arr,log='y',add=T,type='l',col='orange')


# 
# 
# # plot the histogram of empirical meeting times
# mt_max <- max(mt_arr)
# count_arr <- matrix(NA,length(N_arr),mt_max)
# for(n_i in 1:length(N_arr)){
#   for(i in 1:mt_max){
#     count_arr[n_i,i] <- sum(mt_arr[n_i,]==i)
#   }
# }
# prob_arr <- count_arr/dim(mt_arr)[2]
# sd_err <- sqrt(prob_arr*(1-prob_arr)/R)
# prob_arr_u <- prob_arr + 2*sd_err
# prob_arr_l <- prob_arr - 2*sd_err
# 
# plot_df <- data.frame(t(prob_arr))
# plot_df_u <- data.frame(t(prob_arr_u))
# plot_df_l <- data.frame(t(prob_arr_l))
# 
# names(plot_df) <- N_arr
# names(plot_df_u) <- N_arr
# names(plot_df_l) <- N_arr
# 
# plot_df$tau_val <- 1:mt_max
# plot_df_u$tau_val <- 1:mt_max
# plot_df_l$tau_val <- 1:mt_max
# 
# plot_df_m <- reshape2::melt(plot_df,id='tau_val')
# plot_df_m_u <- reshape2::melt(plot_df_u,id='tau_val')
# plot_df_m_l <- reshape2::melt(plot_df_l,id='tau_val')
# 
# 
# names(plot_df_m)[2] <- 'N'
# plot_df <- cbind(plot_df_m,plot_df_m_u$value,plot_df_m_l$value)
# names(plot_df)[4:5] <- c('u','l')
# 
# plot_df_exact <- data.frame(t(p_k_mat))
# names(plot_df_exact) <- N_arr
# plot_df_exact$tau_val <- 1:mt_max
# 
# plot_df_exact_m <- reshape2::melt(plot_df_exact,id='tau_val')
# names(plot_df_exact_m)[2] <- 'N'
# plot_df_exact_m$N <- as.factor(plot_df_exact_m$N)
# 
# g<-ggplot() +
#   geom_bar(data = plot_df_exact_m,aes(x=tau_val,y=value,fill=N),stat='identity',position='dodge')+
#   geom_errorbar(data=plot_df, aes(x=tau_val,y=value,fill=N,ymin=l,ymax=u),position='dodge') +
#   xlab(TeX('$n$'))+
#   scale_x_continuous(breaks=1:100,limits=c(0,10)) +
#   ylab(TeX('$P\\[\\tau = n\\]$'))
# 
# g
# 
# name <- 'lgssm_estimated_mt.pdf'
# ggsave(name,plot=g,width=7,height=5)
# 
# 
# 
# ### compare quantiles
# 
# # compare actual and predicted meeting times
# # for predicted values we need to do a series of integrations over the range of tau appearing in the histograms
# # define the conditional geometric probability
# p_z <- function(z){return(pnorm(-z/sqrt(s2)) + exp(-z+s2/2)*pnorm((z-s2)/sqrt(s2)))}
# # define the integral
# I_z <- function(z){P_z <- (1-p_z(z))^(k-1) * dnorm(z,0,sd=sqrt(s2))}
# # estimate the integral on a grid of values
# k_max <- max(mt_arr)
# k_plt_arr = 1:k_max
# 
# p_k_quantile_mat <- matrix(NA,length(N_arr),k_max)
# pm_quantile_mat <- matrix(NA,length(N_arr),k_max)
# pm_quantile_sd <- matrix(NA,length(N_arr),k_max)
# # ok for each estimate of the variance of the log-likelihood we need to estimate the respective tau probabilities
# for(N_i in 1:length(N_arr)){
#   nparticles <- N_arr[N_i]
#   s2 <- varlogz_arr[N_i]
#   
#   lower = -10*sqrt(s2)
#   upper = 10*sqrt(s2)
#   
#   for(i in 1:length(k_plt_arr)){
#     k <- k_plt_arr[i]
#     p_k_quantile_mat[N_i,i] <- integrate(I_z, lower, upper)$value
#     pm_quantile_mat[N_i,i] <- sum(pmmh_mts[N_i,]>=k)/length(pmmh_mts[N_i,])
#     pm_quantile_sd[N_i,i] <- sqrt(pm_quantile_mat[N_i,i]*(1-pm_quantile_mat[N_i,i]))/sqrt(length(pmmh_mts[N_i,]))
#   }
# }
# 
# 
# delta_offset <- 0.1
# n_data <- 15
# x1 <- c(1:n_data)-delta_offset
# x2 <- c(1:n_data)+delta_offset
# pm_quantile_u <- t(pm_quantile_mat+2*pm_quantile_sd)[1:n_data,c(2,4,6)]
# pm_quantile_l <- t(pm_quantile_mat-2*pm_quantile_sd)[1:n_data,c(2,4,6)]
# 
# 
# ggplot_df1 <- data.frame(t(p_k_quantile_mat[c(2,4,6),1:n_data]),x1)
# names(ggplot_df1)[1:3] <- c((N_arr[2]),(N_arr[4]),(N_arr[6]))
# ggplot_df1_m <- melt(ggplot_df1,id.vars='x1')
# 
# erbar_df1 <- data.frame(t(pm_quantile_mat)[1:n_data,2],pm_quantile_l[1:n_data,1],pm_quantile_u[1:n_data,1],x2)
# erbar_df2 <- data.frame(t(pm_quantile_mat)[1:n_data,4],pm_quantile_l[1:n_data,2],pm_quantile_u[1:n_data,2],x2)
# erbar_df3 <- data.frame(t(pm_quantile_mat)[1:n_data,6],pm_quantile_l[1:n_data,3],pm_quantile_u[1:n_data,3],x2)
# names(erbar_df1)[1:3] <- c('q','q_l','q_u')
# names(erbar_df2)[1:3] <- c('q','q_l','q_u')
# names(erbar_df3)[1:3] <- c('q','q_l','q_u')
# ggplot_df2 <- data.frame(t(pm_quantile_mat[c(2,4,6),1:n_data]),x2)
# ggplot_df2_m <- melt(ggplot_df2,id.vars='x2')
# 
# 
# 
# breaks <- 10^(-10:10)
# minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))
# 
# g <- ggplot() + 
#   geom_point(data=ggplot_df1_m,aes(x=x1,y=value,shape=variable),size=2) +
#   geom_point(data=ggplot_df2_m,aes(x=x2,y=value),size=0.4) +
#   geom_errorbar(data=erbar_df1,aes(x=x2,ymin=q_l,ymax=q_u),width=0.2) +
#   geom_errorbar(data=erbar_df2,aes(x=x2,ymin=q_l,ymax=q_u),width=0.2) +
#   geom_errorbar(data=erbar_df3,aes(x=x2,ymin=q_l,ymax=q_u),width=0.2) +
#   scale_x_continuous(breaks = 0:100,minor_breaks = NULL)+
#   scale_y_log10(breaks =breaks ,minor_breaks=minor_breaks) + 
#   scale_shape_manual(values=c(17, 16, 3))+
#   ylab(TeX('$P\\[\\tau\\geq n\\]$')) +
#   xlab('n') +
#   guides(shape=guide_legend(title="N"))
# g
# name <- ('tailproblgssm.pdf')
# ggsave(name,plot=g,width=7,height=5)
# 
# 


