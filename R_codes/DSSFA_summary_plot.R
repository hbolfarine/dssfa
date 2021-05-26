library(Rcpp)
library(ggplot2)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(latex2exp)

rcpp_inc_plot<- '
using namespace Rcpp;
using namespace arma;
# '

src.fanc.plot.cov = '
cube Omegacov   = as<cube>(Mat1);
cube Omegakmax  = as<cube>(Mat2);

int M = size(Omegacov)(2);
int kmax = size(Omegakmax)(2);
arma::mat Omega_lamb(M,kmax);
double val;
double sign;

for(int i = 0; i < kmax; i++){
  log_det(val,sign,Omegakmax.slice(i));
  for(int j = 0; j < M; j++){
    Omega_lamb(j,i) = val+trace(inv_sympd(Omegakmax.slice(i))*(Omegacov.slice(j)));
  }
}

return(wrap(Omega_lamb));
'
fitcpp.post.plot.cov = cxxfunction(signature(Mat1 = "numeric",Mat2 = "numeric"), src.fanc.plot.cov, plugin='RcppArmadillo', rcpp_inc_plot)


# rcpp_inc_plot<- '
# using namespace Rcpp;
# using namespace arma;
# # '

src.fanc.plot = '
cube load      = as<cube>(Mat1);
mat  uniq      = as<mat>(Mat2);
cube Omegakmax  = as<cube>(Mat3);

int M = size(load)(2);
int kmax = size(Omegakmax)(2);
arma::mat Omega_lamb(M,kmax);
double val;
double sign;

for(int i = 0; i < kmax; i++){
  log_det(val,sign,Omegakmax.slice(i));
  for(int j = 0; j < M; j++){
    Omega_lamb(j,i) = val+trace(inv_sympd(Omegakmax.slice(i))*(load.slice(j)*load.slice(j).t()+diagmat(uniq.row(j))));
  }
}

return(wrap(Omega_lamb));
'
fitcpp.post.plot = cxxfunction(signature(Mat1 = "numeric",Mat2 = "numeric",Mat3 = "numeric"), src.fanc.plot, plugin='RcppArmadillo', rcpp_inc_plot)


summary.plot.dssfa = function(list.bfa,CI = NULL, cov.dssfa = F, reg = T, text_main_prior = "", col.lambda = 1, legend.pos = "right"){
  
  if(is.null(CI)){
    ic.1 = 0.025
    ic.2 = 0.975
    CI.text = "95"
  }else{
    ic.1 = CI[[2]][1]
    ic.2 = CI[[2]][2]
    CI.text = CI[[1]]
  }
  
  # Quantis
  q1 = quantile(list.bfa$fit.post.fanc,probs = c(ic.1,ic.2))[1]
  q2 = quantile(list.bfa$fit.post.fanc,probs = c(ic.1,ic.2))[2]
  
  # Select number of factors and lambda
  k.sel = min(which(lapply(list.bfa$fit.sel, function(x) length(which(x < q2))) != 0))
  lamb.sel = lapply(list.bfa$fit.sel, function(x) list.bfa$length.lambda-length(which(x < q2))+1)
  lamb.sel = lamb.sel[[k.sel]]
  
  fit.list.plot = list()
  for(i in 1:list.bfa$k.max){
    if(length(dim(list.bfa$list.sel[[i]])) == 3){
      dim.size = dim(list.bfa$list.sel[[i]])[3]
      Omega.temp = array(0,dim = c(list.bfa$p,list.bfa$p,dim.size))
      for(j in 1:dim.size){
        beta.max = as.matrix(list.bfa$list.sel[[i]][,,j])
        uniq.sel = list.bfa$uniq.full[[i]][j,]
        Omega.temp[,,j] = beta.max%*%t(beta.max)+diag(uniq.sel)
      }
    }else if(length(dim(list.bfa$list.sel[[i]])) != 3){
      dim.size = 1
      Omega.temp = array(0,dim = c(list.bfa$p,list.bfa$p,1))
      beta.max = as.matrix(list.bfa$list.sel[[i]])
      uniq.sel = list.bfa$uniq.full[[i]]
      Omega.temp[,,1] = beta.max%*%t(beta.max)+diag(uniq.sel)
    }
    if(cov.dssfa == T){
      fit.list.plot[[i]] = as.matrix(fitcpp.post.plot.cov(list.bfa$post.cov,Omega.temp))
    }else{
      fit.list.plot[[i]] = as.matrix(fitcpp.post.plot(list.bfa$post.loadings,list.bfa$Sigma.post,Omega.temp))
    }
    colnames(fit.list.plot[[i]]) = ((list.bfa$length.lambda-dim.size)+1):list.bfa$length.lambda
  }
  
  ##################################
  
  temp.melt = lapply(fit.list.plot, melt)
  
  for(i in 1:list.bfa$k.max){
    temp.melt[[i]]$kmax = as.numeric(i)
  }
  
  temp = bind_rows(temp.melt) #dplyr
  temp = data.frame(temp)
  
  max.lamb = max(temp$Var2)
  temp$Var2 = factor(temp$Var2,levels = 1:max.lamb,ordered = T)
  
  temp.melt2 = lapply(list.bfa$fit.sel,melt)
  for(i in 1:list.bfa$k.max){
    temp.melt2[[i]]$kmax = as.numeric(i)
    dim.size = dim(temp.melt2[[i]])[1]
    temp.melt2[[i]]$lamb = ((list.bfa$length.lambda-dim.size)+1):list.bfa$length.lambda
  }
  
  temp2 = bind_rows(temp.melt2)
  # temp.melt2$lamb = factor(rev(temp.melt2[[i]]$lamb), levels = 1:max.lamb,ordered = T)
  
  mygreys = colorRampPalette(brewer.pal(9,"Greys"))(max.lamb + 5)
  mygreys = mygreys[-(1:5)]

  min = min(temp2$lamb)
  g = ggplot() + geom_point(data = temp2, aes(factor(kmax),value, col = factor(lamb, labels = (max.lamb-min):0)), size = 1) + theme_bw() + theme(legend.position = "none") + labs(col=TeX("$\\lambda_{i}$")) +
    scale_color_manual(values = mygreys, aesthetics = c("col")) + theme(legend.position = legend.pos) +
    guides(colour = guide_legend(ncol = col.lambda))
  
  temp.data.reg = temp[which(temp$kmax == k.sel & temp$Var2 == lamb.sel),]
  temp.point.reg = temp2[which(temp2$kmax == k.sel & temp2$lamb == lamb.sel),]
  
  for(i in 1:max.lamb){
    
    temp.data = temp[which(temp$Var2 == i),]
    temp.point = temp2[which(temp2$lamb ==  i),]
    g = g + geom_violin(data  = temp.data,aes(factor(kmax), value), alpha = 0, size = 0.4, col = mygreys[i]) + geom_point(data = temp.point, aes(kmax,value,col = factor(lamb)), col = mygreys[i])
    
    if(reg == T & i == lamb.sel){
      g = g + geom_violin(data = temp.data.reg, aes(factor(kmax),value),col = "red", alpha = 0, size = 0.7) + geom_point(data = temp.point.reg, aes(kmax,value), col = "red")
    }    
    
  }
  
  if(reg == T){
    lamb.sel = max.lamb-lamb.sel
    text.lambda = paste("$\\lambda_{",lamb.sel,"}$", sep = "")
    text.main = paste("DSSFA Model Summary,",text_main_prior," $k_{dssfa} = $",list.bfa$k.max,", ",CI.text,"% CI, selected (",text.lambda,", $k = $",k.sel,")",sep = "")
  }else{
    text.main = paste("DSSFA Model Summary,",text_main_prior," $k_{dssfa} = $",list.bfa$k.max,", ",CI.text,"% CI",sep = "")  
  }
  g = g + geom_hline(yintercept = q2, lty = "dashed", col = "darkgray") + geom_hline(yintercept = q1, lty = "dashed", col = "darkgray") + 
    labs(title = TeX(text.main), x = "Number of Factors, k", y = TeX("loss function")) #+ coord_flip()
      # labs(title = TeX("DSSFA Model Summary, normal prior, $k_{dssfa}=24$, 95% CI"), x = "Number of Factors, k", y = TeX("$L_{\\lambda}^{k}(\\Omega$ | $\\Omega_{\\lambda}^{k})$")) #+ coord_flip()
  g
  
}

