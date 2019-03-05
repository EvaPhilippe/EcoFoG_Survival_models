

functions{
  real log_F(
    int r, //recrutement
    int k, //date pour laquelle on calcule 
    int[] dates, //vecteur des dates d'inventaires, en commençant à 1 et pas 1992
    row_vector beta, //vecteur des paramètres
    vector[] x){ //vecteur des covariables par année pour un individu
    real ret ;
    real ret2 ;
    ret=0 ;
    for (j in dates[r]:(dates[k-1]-1)){
      ret-=exp(beta*x[j]) ;
    }
    ret2=0 ;
    for (j in dates[k-1]:(dates[k]-1)){
      ret2-=exp(beta*x[j]) ;
    }
    ret+=log_diff_exp(0, ret2) ;
    return ret ;
  }
  real log_likelihood(
    int r, //recrutement
    int m, //mort 
    int[] dates, //vecteur des dates d'inventaire
    row_vector beta, //vecteur des paramètres
    vector[] x, //vecteur des covariables par année
    int cens //indicatrice de censure
    ){
      if (cens==0){
        return log_F(r, m, dates, beta, x);
      }
      else{
        real ret ;
        ret=1 ;
        for (k in (r+1):m){
          ret-=exp(log_F(r, k, dates, beta, x)) ;
        }
      return log(ret) ;
      }
    }
}


data {
  int<lower=1> D ; // nombre d'inventaires
  int Dates[D] ; // Dates[j] correspond à la date de l'inventaire j, en commençant à 1 et pas 1992
  int<lower=1> I ; //nombre d'individus
  int<lower=0> P; //nombre de covariables+1
  int<lower=1> L; //nombre total d'années sur la période étudiée
  vector[P] X[I,L] ; // array des covariables par individu et par année, X[i,j, c] correspond à la covariable c de l'individu i à l'inventaire j
  int Recr[I] ; //Recr[i] correspond à l'indice de l'inventaire de recrutement de l'individu i
  int Mort[I] ; //Mort[i] correspond soit à l'inventaire auquel l'individu i est observé mort, soit au dernier inventaire auquel il a été vu vivant, s'il est censuré à droite
  int Cens[I] ; //Cens[i]=1 si l'individu est censuré (sa mort n'est pas observée), 0 sinon
}


parameters{
  row_vector[P] Beta ; //vecteur des paramètres, Beta[1] est l'intercept
}


model{
  for (i in 1:I){
    target+= log_likelihood(
      Recr[i], 
      Mort[i], 
      Dates,
      Beta,
      X[i],
      Cens[i]) ;
  }
}

//generated quantities{
//  row_vector[P] Lambda ;
//  Lambda=exp(Beta) ; //On veut otenir les estimations de lambda
//}
