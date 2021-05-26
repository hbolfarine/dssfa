############################################################################################
#
# "Decoupling Shrinkage and Selection in Gaussian Linear Factor Analysis", 2020
# - Summary plots, Selected loadings
#
############################################################################################
#
# Author : Henrique Bolfarine
#          Institute of Mathematics and Statistics
#          University of São Paulo
#          Email : bolfarin@ime.usp.br
#
############################################################################################
# 
# Description: Infere the number of factors in the international exchange dataset

#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#

# Function that selects the loadings, given the CI and the number of factors.
sel.model.dssfa = function(list.bfa,CI = NULL, k.sel = NULL, fast = F){
  
  if(length(CI) == 0){
    # 95%  credible interval defaut
    ic.1 = 0.025
    ic.2 = 0.975
  }else{
    ic.1 = CI[[2]][1]
    ic.2 = CI[[2]][2]
  }
  
  if(is.null(k.sel)){
    k.sel = 1:list.bfa$k.max
  }
  
  # Get quantiles from the credible interval.
  q1 = quantile(list.bfa$fit.post.fanc,probs = c(ic.1,ic.2))[1]
  q2 = quantile(list.bfa$fit.post.fanc,probs = c(ic.1,ic.2))[2]
  
  sel.model     = list()
  mat.sel.ksel  = vector("list",length(k.sel))
  uniq.sel.ksel = vector("list",length(k.sel))
  phi.sel.ksel  = vector("list",length(k.sel))
  mat.sel.CI    = vector("list",list.bfa$k.max)
  uniq.sel.CI   = vector("list",list.bfa$k.max)
  phi.sel.CI    = vector("list",list.bfa$k.max)
  
  # mod.sel = rep(0,length(k.sel))
  # k.sel.1 = rep(0,list.bfa$k.max)
  # k.sel.1[k.sel] = k.sel
  
  l = 1 # count the number of select factors
  if(fast){
    for(i in k.sel){
      if(sum(which(list.bfa$fit.sel[[i]] < q2)) == 0){
        # cat("aqui11")
        mat.sel.ksel[[l]]  = NULL
        uniq.sel.ksel[[l]] = NULL
        if(list.bfa$cor.factor){
          phi.sel.ksel[[l]] = NULL
        }
        l = l+1
      }else{
        # cat("aqui12")
        min.sel = min(which(list.bfa$fit.sel[[i]] < q2))
        mat.sel.ksel[[l]]  = list.bfa$betas.full[[i]]
        uniq.sel.ksel[[l]] = list.bfa$list.uniqueness[[i]]
        if(list.bfa$cor.factor){
          phi.sel.ksel[[l]] = list.bfa$phi.full[[i]]
        }
        l = l+1
      }
    }
  }else{
    for(i in k.sel){
      if(sum(which(list.bfa$fit.sel[[i]] < q2)) == 0){
        # cat("aqui11")
        mat.sel.ksel[[l]]  = NULL
        uniq.sel.ksel[[l]] = NULL
        if(list.bfa$cor.factor){
          phi.sel.ksel[[l]] = NULL
        }
        l = l+1
      }else{
        if(length(dim(list.bfa$list.sel[[i]])) != 3){
          # cat("aqui12")
          min.sel = min(which(list.bfa$fit.sel[[i]] < q2))
          mat.sel.ksel[[l]]  = list.bfa$list.sel[[i]]
          uniq.sel.ksel[[l]] = list.bfa$uniq.full[[i]]
          if(list.bfa$cor.factor){
            phi.sel.ksel[[l]] = list.bfa$list.bfa$phi.list[[i]]
          }
          l = l+1
        }else{
          # cat("aqui13")
          min.sel = min(which(list.bfa$fit.sel[[i]] < q2))
          mat.sel.ksel[[l]]  = list.bfa$list.sel[[i]][,,min.sel]
          if(list.bfa$length.lambda < 2){
            uniq.sel.ksel[[l]] = list.bfa$uniq.full[[i]]
          }else{
            uniq.sel.ksel[[l]] = list.bfa$uniq.full[[i]][min.sel,]
          }
          if(list.bfa$cor.factor){
            if(i == 1){
              phi.sel.ksel[[l]] = 1
            }else{
              phi.sel.ksel[[l]] = list.bfa$phi.list[[i]][,,min.sel]
            }
          }
          l = l+1
        }
      }
    }
  }
  #----------------------------------------------------#
  
  sel.model$mat.sel.ksel  = mat.sel.ksel
  sel.model$uniq.sel.ksel = uniq.sel.ksel
  if(list.bfa$cor.factor){
    sel.model$phi.sel.ksel = phi.sel.ksel
  }
  
  return(sel.model)
  
}

make.plot.mat = function(mat,col.names,fact.name){
  rownames(mat) = col.names
  colnames(mat) = fact.name
  mat = apply(mat, 2, rev)
  longmat = melt(mat)
  
  names(longmat)[2] = "factors"
  names(longmat)[1] = "variable"
  
  g1 = ggplot(longmat, aes(x = factors, y = variable)) + geom_tile(aes(fill = value),colour = "grey20")
  
  g2 = g1 + scale_fill_gradient2(low = "#800000", high = "#191970", mid = "white")+ theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    # panel.grid.major = element_blank(), 
    text = element_text(size = 12),
    panel.border = element_blank(),
    panel.background = element_blank(), 
    axis.ticks = element_blank(),
    # axis.text = element_blank(), 
    legend.title = element_text(),
    plot.title = element_text(hjust = 0.5)) + labs(fill = " ")
  g2
}


#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
