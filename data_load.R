#data input
library(dplyr)
library(splines)
#functions--------
#serial interval distribution
weibull_offset=function(shape, scale, d, week=T){
  k=ifelse(week,7,1)
  k0=ifelse(week,1,5)
  lapply(0:(50*k0), function(x){
    pweibull(x*k+k-d,shape, scale=scale)-
      pweibull(x*k-d,shape, scale=scale)
  })%>%unlist()
}

#age group aggregation
aggregateagegroups <- function(p, newageinterval) {
  if (!all(newageinterval %in% p$ageinterval)) stop("new age intervals are not subset of existing intervals")
  ind <- match(newageinterval, p$ageinterval)
  if(is.matrix(p$cases)){
    agg_case=as.matrix(do.call('rbind',lapply(1:nrow(p$cases), function(j){
      agg_case=lapply(2:length(ind), function(k){
        sum(p$cases[j,ind[k-1]:(ind[k]-1)])
      })%>%unlist
      agg_case=c(agg_case,sum(p$cases[j,-c(1:ind[length(ind)]-1)]) )
    })))
  }else{
    agg_case=lapply(2:length(ind), function(k){
      sum(p$cases[ind[k-1]:(ind[k]-1)])
    })%>%unlist
    agg_case=c(agg_case,sum(p$cases[-c(1:ind[length(ind)]-1)]) )
  }
  list(cases = agg_case, casecategories = p$casecategories, 
       ageinterval = newageinterval, misc = p$misc)
}

#defining ContactMatrix structure
ContactMatrix <- function(matrix, ageinterval, susceptibility, parameters, misc) {
  list(matrix = matrix, ageinterval = ageinterval, susceptibility = susceptibility, 
       parameters = parameters, misc = misc)
}

#contact matrix derivation
newcontact <- function(surveyname, ageinterval, year=2024,filter = NULL,latest_pop) {
  pop=data.frame(lower.age.limit=ageinterval, 
                 population=latest_pop)
  smixr <- socialmixr::contact_matrix(surveyname, survey.pop=pop,countries = 'Zimbabwe', 
                                      age.limits = ageinterval, filter = filter, symmetric = TRUE, 
                                      return.demography = TRUE, estimated.contact.age = "sample",weigh.age = T)
  cmt <- smixr$matrix/2
  ContactMatrix(cmt, ageinterval, rep(1,length(ageinterval)), list(),list(issynthetic = F)) 
}

#adjust contact matrix based on updated population
contact_adjust=function(cm,new_pop, old_pop, scale=T){
  if(scale){
    c0=sum(old_pop)/sum(new_pop)
  }else{
    c0=1
  }
  lapply(1:nrow(cm), function(i){
    (cm[i,]*new_pop)/old_pop*c0
  })%>%do.call('rbind',.)%>%as.matrix
}

#aggregating zimbabwe population based on given age interval
zwe_pop=read.csv('data/ZWE_pop_2012.csv')
zwe_pop_combine=function(ageint){
  unlist(lapply(seq_along(ageint), function(k){
    min_age=ifelse(k>1,ageint[k-1],-1)
    max_age=ageint[k]
    idx=which(zwe_pop$age<=max_age&zwe_pop$age>min_age)
    sum(zwe_pop$pop[idx])
  }))
}

#Brunudi data-------
#contact & population
new_pop=read.csv('data/BDI_pop_2024.csv') #for each age
ageint=c(0:3*5,2:5*10)
bdi_pop=lapply(seq_along(ageint), function(k){
  idx1=which(new_pop$age_lower==ageint[k])
  if(k<length(ageint)){
    idx2=which(new_pop$age_lower==ageint[k+1])-1
    idx=idx1:idx2
  }else{
    idx=idx1:nrow(new_pop)
  }
  colSums(new_pop[idx,-1])/1e3
})%>%do.call('rbind',.)%>%as.matrix
zwe_pop1=zwe_pop_combine(ageint)/1e3
zwe2012_home <- newcontact(zimbabwe_survey, ageint, "ZWE", 2012,filters[[3]], zwe_pop1)$matrix
bdi_cm=contact_adjust(zwe2012_home,rowSums(bdi_pop),zwe_pop1, scale=T)

#serial interval
w=weibull_offset(1.63,18.29,1)

#sexually active group
sex_idx=which(ageint%in%c(3:9*5))

n.age=length(ageint)
N0=3 #days without observations
n.t=nrow(case_mat)+N0
trange=1:n.t
df=6

#data for stan model
data_list=list(
  N0=N0, N=n.t,
  nage=n.age,
  n_sexage=length(sex_idx), sex_idx=sex_idx,
  cases=case_mat,
  cm=bdi_cm,
  v_cover=c(rep(0,n.age-2),0.5,1),
  n_SI=length(w), SI=w,
  na=bdi_pop,
  sex_p=rep(0.1,2),
  K=df, 
  B=bs(1:n.t, df=df,degree=3, intercept=T),
  temp=1/log((n.t-N0)*n.age*2-2*N1)
)



