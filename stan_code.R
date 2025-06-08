#estimating using case pyramids only
stan_code <- "
data {
  int<lower=1> N; // days of transmission
  int<lower=0> N0; // number of days without case pyramids
  int<lower=1> nage; // number of age categories
  array[N-N0, nage*2] int cases; // reported case counts by age and sex
  int n_sexage; // number of sexually active age groups
  array[n_sexage] int sex_idx; // index of sexually active age groups
  matrix[nage, nage] cm; //contact matrix
  int n_SI;
  vector[n_SI] SI; // fixed pre-calculated SI 
  vector[nage] v_cover; // vaccine coverage
  vector[2] sex_p; //overall highly sexually active proportion
  matrix[nage,2] na; //population (thousand)
  int<lower=1> K;       // number of spline basis functions
  matrix[N, K] B;       // spline basis matrix
  real<lower=0> temp; //temperature 1/log(#observations)
}

parameters {
  vector<lower=1e-6>[nage*4] y0;
  matrix<lower=1>[n_sexage,2] paqa_raw; //sexually active proportion para
  vector[K] beta_0;
  matrix[K,2] beta_s;
  matrix[K,2] beta_c;
  vector<lower=0>[2] sigma_s;
  vector<lower=0>[2] sigma_c;
  vector[2] mu_s;
  vector[2] mu_c;
  real<lower=0> ve_log;
}

transformed parameters {
    real ve;
    vector[nage] sus;
    vector[nage*4] sus1;
    vector<lower=0>[N] c;
    matrix<lower=0>[N,2] sus_c; //children susceptibility (0-4, 5-9)
    matrix<lower=0>[N,2] vw; //w_F, w_M
    matrix<lower=0>[N, nage*4] inf=rep_matrix(1e-6,N,nage*4); //total infectiousness
    matrix<lower=0>[nage*4,nage*4] sus_m;
    matrix<lower=0>[nage*4,nage*4] cmt=rep_matrix(0,nage*4,nage*4); //transpose of the augmented contact matrix
    matrix<lower=0>[nage*4,nage*4] ngm=rep_matrix(0,nage*4,nage*4); //next generation matrix
    matrix<lower=0>[nage,2] paqa=rep_matrix(0,nage,2); //sexually active proportion
    matrix<lower=0>[nage,2] napaqa=rep_matrix(0,nage,2); //sexually active count
    matrix<lower=0>[nage,2] Q_vec=rep_matrix(1, nage, 2); //for sexual mat construction
    matrix<lower=0>[nage, nage] sex_mat=rep_matrix(0,nage, nage); //temporary sexual matrix
    vector<lower=0>[nage] sex_vec=rep_vector(0,nage); //temporary sexual vector
    vector<lower=0>[2] sum_napaqa=rep_vector(0,2); //sum of sexually active pop
    vector<lower=0>[nage] p_male=rep_vector(0,nage); //sum of sexually active pop
    vector<lower=0>[2] vw0;
    matrix<lower=0>[N,nage*4] prediction;
    matrix<lower=0>[N,nage*2] I;
    vector<lower=0>[N] I_total; //total case counts
    p_male=na[,1]./(na[,1]+na[,2]);
    //construct sexual transmission
    for(i in 1:2){
      paqa[sex_idx,i]=pow(1 - paqa_raw[1:n_sexage,i],2)/2;
      napaqa[1:nage,i]=paqa[1:nage,i].*na[1:nage,i];
      sum_napaqa[i]=sum(napaqa[1:nage,i]);
    }
  //initial case counts
  prediction[1,]=to_row_vector(y0);
  // expected case counts on day t
  c=exp(B*beta_0);
  for(i in 1:2){
    sus_c[,i]=to_vector(exp(B*beta_c[,i]));
    vw[,i]=to_vector(exp(B*beta_s[,i]));
  }
  for (i in 2:N) {
    vw0[1]=vw[i,1] / 2.175; //coeffient for Sigma_MF
    vw0[2]=vw0[1] * (sum_napaqa[2] / sum_napaqa[1]); //coeffient for Sigma_FM
    for(m in 1:2){ //two sexes
      Q_vec[1:nage,2]=paqa[1:nage,3-m];
      for(j in 1:2){ //two sexually active levels
        sex_mat=rep_matrix(0,nage, nage);
        sex_vec=vw[i,m]/sum_napaqa[m]*(na[1:nage,m].*paqa[1:nage,m]);
        for(k in sex_idx){ //sexually active age groups
          sex_mat[1:nage,k]=sex_vec*Q_vec[k,j];
        }
        for(k1 in 1:nage){
          for(k2 in 1:nage){
            cmt[(m-1)*2*nage+k1,(2-m)*nage*2+(j-1)*nage+k2]=sex_mat[k1,k2];
          }
        }
       }
      //community contact
      for(j in 1:2){
        sex_mat=cm'*diag_matrix(p_male);
        for(k1 in 1:nage){
          for(k2 in 1:nage){
            cmt[(m*2-1)*nage+k1,(j-1)*nage+k2]=sex_mat[k1,k2];
          }
        }
        sex_mat=cm'*diag_matrix(1-p_male);
        for(k1 in 1:nage){
          for(k2 in 1:nage){
            cmt[(m*2-1)*nage+k1,nage*2+(j-1)*nage+k2]=sex_mat[k1,k2];
          }
        }
      }
    }
    ve=exp(-ve_log); //vaccine effectiveness
    sus=1-v_cover*ve;
    sus[1:2]=to_vector(sus_c[i,1:2]).*sus[1:2];
    for(j in 0:3){
      sus1[(1+nage*j):(nage+nage*j)]=sus;
    }
    sus_m=diag_matrix(sus1);
    ngm=sus_m*cmt;
    for(j in 1:(4*nage)){
      inf[i,j]=dot_product(reverse(SI[1:(i-1)]), prediction[1:(i-1),j]);
    }
    ngm=ngm*c[i];
    prediction[i,] = to_row_vector(ngm * inf[i,]');
  }
  for(i in 1:N){
    I[i,1:nage]=prediction[i,1:nage]+prediction[i,(1+nage):(2*nage)];
    I[i,(1+nage):(nage*2)]=prediction[i,(1+2*nage):(3*nage)]+prediction[i,(1+3*nage):(4*nage)];
  }
  //total case counts
  for(i in 1:N){
    I_total[i]=sum(I[i,]);
  }
}
model {
  ve ~ normal (0.75,0.15);
  y0 ~ normal(0,10);
  sigma_c~normal(0,10);
  sigma_s~normal(0,10);
  mu_c~normal(0,10);
  mu_s~normal(0,10);
  c ~ normal(0.15,0.2);
  for(i in 1:2){
    beta_c[,i]~normal(mu_c[i],sigma_c[i]);
    beta_s[,i]~normal(mu_s[i],sigma_s[i]);
  }
  //for those with total case counts
  for(i in 1:2){
    sum_napaqa[i]~normal(sex_p[i]*sum(na[sex_idx,i]),1);
  }
  for(i in (N0+1):N){
    cases[i-N0,] ~ poisson(I[i,]);
  }
}
"
