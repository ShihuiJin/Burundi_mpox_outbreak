#functions---------
estim=function(dat){
  if(is.vector(dat)){
    c(mean(dat),quantile(dat,c(0.5,0.025,0.975,0.25,0.75)))
  }else{
    as.matrix(do.call('rbind',lapply(1:ncol(dat), function(i){
      c(mean(dat[,i]),quantile(dat[,i],c(0.5,0.025,0.975,0.25,0.75)))
    })))
  }
}

construct_sex=function(na,pa, coef, Q_vec=NULL){
  vec=coef*na*pa/sum(na*pa)
  if(is.null(Q_vec)) Q_vec=rep(1, length(na))
  as.matrix(do.call('cbind',lapply(1:length(pa), function(i){
    if(pa[i]==0){
      rep(0,length(pa))
    }else{
      vec*Q_vec[i]
    }
  })))
}

ngm_Rt=function(sus_c, R0_scale, cover, cm=bdi_cm, sex=T, 
                ve=0.8, paqa=NULL,vw=NULL,na=bdi_pop, fit_eigen=T){
  suscept=1-cover*ve
  suscept[1:2]=sus_c
  if(sex==F){ #non-sex model: clade 1a?
    ngm=diag(suscept)%*%t(cm)
    fit=eigen(ngm)
    vec=abs(fit$vectors[,1])
    vec=vec/sum(vec)
    val=max(abs(fit$values))
    list(val=val*R0_scale,vec=vec,ngm=ngm)
  }else{ #with sexual transmission
    sum_napaqa <- sapply(seq_along(paqa), function(x) sum(na[,x] * paqa[[x]]))
    v0 <- vw[1] / ((sqrt(12^2 + 22^2) / (17.61 + 5.51))^2 + 1) #coeffient for Sigma_MF
    w0 <- v0 * (sum_napaqa[2] / sum_napaqa[1]) #coeffient for Sigma_FM
    s_cmt=list(
      construct_sex(na[,1],paqa[[1]],vw[1]), #S_MF
      construct_sex(na[,2],paqa[[2]],vw[2]), #S_FM
      construct_sex(na[,1],paqa[[1]],v0, paqa[[2]]), #Sigma_MF*Q
      construct_sex(na[,2],paqa[[2]],w0, paqa[[1]]) #Sigma_FM*P
    )
    n.col=length(paqa[[1]])
    null_mat=matrix(0,n.col,n.col)
    p_male=na[,1]/rowSums(na)
    contact_mats=as.matrix(do.call('cbind',lapply(1:4, function(i){
      if(i<=2){
        t(cm)%*%diag(p_male)
      }else{
        t(cm)%*%diag(1-p_male)
      }
    })))
    sex_mat=rbind(
      cbind(null_mat,null_mat, s_cmt[[1]],s_cmt[[3]]),
      contact_mats,
      cbind(s_cmt[[2]],s_cmt[[4]],null_mat,null_mat),
      contact_mats
    )%>%as.matrix%>%t()
    suscept1=rep(suscept, 4)
    ngm=diag(suscept1)%*%t(sex_mat)
    if(fit_eigen){
      fit=eigen(ngm)
      vec0=abs(fit$vectors[,1]) #unscaled
      vec=vec0/sum(vec0)
      vec1=vec[c(1:n.age, 1:n.age+n.age*2)]+vec[n.age+c(1:n.age, 1:n.age+n.age*2)]
      val=max(abs(fit$values))
      list(val=val*R0_scale, vec_combined=vec1, vec=vec0, ngm=ngm)
    }else{
      ngm
    }
  }
}
#analysis--------
temp=fit$draws(format = 'matrix')
trange=(N0+1):n.t
#VE
estim(temp[,'ve'])

#sexually active %
lapply(1:2, function(i){
  estim(temp[,paste0('paqa[',sex_idx,',',i,']')]*100)%>%round(2)
})

#w_F & w_M
para_vw=lapply(1:2, function(i){
  estim(temp[,paste0('vw[',trange,',',i,']')])
})

#relative susceptibility
para_sus=lapply(1:2, function(t){
  estim(temp[,paste0('sus_c[',trange,',',t,']')])
})

#normalizing constant
para_c=estim(temp[,paste0('c[',trange,']')])

#Rt
eign_estim=lapply(trange, function(t){
  sus_c=temp[,paste0('sus_c[',t,',',1:2,']')]
  R0_scale=temp[,paste0('c[',t,']')]
  vw_t=temp[,paste0('vw[',t,',',1:2,']')]
  paqa_list=lapply(1:nrow(temp), function(i){
    list(
      temp[i,paste0('paqa[',1:n.age,',1]')]%>%as.vector(),
      temp[i,paste0('paqa[',1:n.age,',2]')]%>%as.vector()
    )
  })
  if(t>1){
    inf_t=lapply(1:(4*n.age), function(j){
      temp[,paste0('prediction[',1:(t-1),',',j,']')]%*%w[(t-1):1]
    })%>%do.call('cbind',.)%>%as.matrix()
  }
  mclapply(1:nrow(temp), function(i){
    ngmfit=ngm_Rt(sus_c[i], R0_scale[i], v_cover, bdi_cm, paqa=paqa_list[[i]],vw=vw_t[i,])
    ngm=ngmfit$ngm
    if(t>1){
      child_parent=lapply(1:2, function(k){
        newinf=(ngm[k,]+ngm[k+n.age,]+ngm[k+n.age*2,]+ngm[k+n.age*3,])*inf_t[i,]
        sum(newinf[sex_idx_all])/sum(newinf)
      })%>%unlist
    }else{
      child_parent=rep(0,2)
    }
    c(ngmfit$val,child_parent)
  },mc.cores = detectCores())%>%do.call('rbind',.)%>%as.matrix()
})
Rt=lapply(eign_estim, function(x) estim(x))%>%do.call('rbind',.)%>%as.matrix()
#% sexual transmission
p_sex=lapply(trange, function(t){
  x1=temp[,paste0('prediction[',t,',',c(1:n.age, 1:n.age+n.age*2),']')]%>%rowSums()
  x2=temp[,paste0('I_total[',t,']')]
  x1/x2*100
})%>%do.call('cbind',.)%>%as.matrix()%>%estim
#% infections aged 0-4 or 5-9 attributable to infectors aged 15-49
child_parent=lapply(1:2+2, function(i){
  lapply(eign_estim, function(x) estim(x)[i,]*100)%>%do.call('rbind',.)%>%as.matrix()
})



