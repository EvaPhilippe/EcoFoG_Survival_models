
## fun_data --------------

#Première fonction utilisée avec le modèle 1 (données uniquement pour les dates d'inventaires) 
#pour passer de X un array à trois dimensions (nb individus, nb d'inventaires, nb de covariables testées+1) 
#à la liste data

fun_data=function(X){
  return(list(
    D=D,
    Dates=Dates_inv,
    I=I,
    P=dim(X)[3],
    X=X,
    Recr= data_indiv$index_recr,
    Mort= data_indiv$index_mort,
    Cens= data_indiv$censure
  ))
}

## fun_data2 -------------

# Fonction pour le Modèle 2, pour le vecteur X des covariables il y a cette fois 
#un vecteur de covariables par individu et par année entre 1992 et 2016 
#(pas seulement pour les années d'inventaires).

#Fonction pour produire la liste data à partir d'un tableau de données X et d'une date (d'inventaire) 
#de départ à partir de laquelle on commence à prendre en compte les inventaires et les mesures.


fun_data2=function(X, start=1992){
  Dates<-dates$annee
  Dates_start<-Dates[which(Dates>=start)]
  data_indiv_start<-data_indiv %>% 
    filter(mort>start) %>% 
    mutate(recrutement=sapply(recrutement, function(x) max(x, start))) %>% 
    mutate(index_recr=sapply(recrutement, function(x) which(Dates_start==x))) %>% 
    mutate(index_mort=sapply(mort, function(x) which(Dates_start==x)))
  return(list(
    D=length(Dates_start),
    Dates=Dates_start-(start-1),
    I=nrow(data_indiv_start),
    P=dim(X)[3],
    L=dim(X)[2],
    X=X,
    Recr= data_indiv_start$index_recr,
    Mort= data_indiv_start$index_mort,
    Cens= data_indiv_start$censure
  ))
}

## fun_data_fill -----------------

#Fonction pour remplir les données (par interpolation linéaires) pour une covariable 
#dont les données n'ont été mesurées que à certaines dates, 
#l'année d'origine est toujours la première année où on a des données 
#(que l'on suppose coïncider avec une année d'inventaire) 
#et on suppose qu'on a toujours des données jusqu'à la fin de l'étude (2016).
#Renvoie un tableau au format long


fun_data_fill <- function(
  tabl_data #tableau des données à utiliser, au format long avec comme colonnes ID, annee (valeurs possibles=celles dans dates_data+première date-1), data. On suppose que tous les ID et toutes les années sont remplies.
){
  dates_data<-sort(as.numeric(unique(tabl_data$annee)))#vecteur des dates auxquelles on a pris des mesures
  D<-length(dates_data)
  I<- length(unique(tabl_data$ID))
  L<- dates_data[D]-dates_data[1]+1 #nombre total de dates à remplir 
  if (nrow(tabl_data) != D*I) stop("There should be ", D, " observations for each of the ", I, " individuals")
  ret<- data.frame(ID=rep(unique(tabl_data$ID), each=L), annee=rep(dates_data[1]:dates_data[D],I)) %>% 
    full_join(tabl_data, by=c("ID", "annee"))
  for (i in 1:I){
    for (k in 1:(D-1)){
      l<-dates_data[k]-dates_data[1]+1 #(pour avoir des indices et pas les années)
      r<-dates_data[k+1]-dates_data[1]+1
      for (j in ((l+1):(r-1))){
        ret$data[L*(i-1)+j]<-ret$data[L*(i-1)+l]+(j-l)/(r-l)*(ret$data[L*(i-1)+r]-ret$data[L*(i-1)+l])
      }
    }
  }
  
  return (ret)
}


## fun_scale ----------------

#Fonction pour centrer-réduire les tableaux de données 
#(en prenant la moyenne et l'écart-type sur les moyennes par individu), 
#de format long à format long.


fun_scale<-function(tabl_data){ #tableau des données à utiliser, au format long avec comme colonnes ID, annee, data.
  moy_indiv<-tabl_data %>%
  group_by(ID) %>%
  summarize(data_moy=mean(data, na.rm=TRUE))
m<-mean(moy_indiv$data_moy, na.rm=TRUE) #moyenne
sd<-sd(moy_indiv$data_moy, na.rm=TRUE) #écart-type
return(tabl_data %>% 
         mutate(data=(data-m)/sd))
}

## fun_data3 --------------------

#Fonction encore plus générique que fun_data2, 
#qui renvoie une liste à utiliser comme data d'entrée pour le modèle script_gen2, 
#à partir d'une liste de tableaux de données centrées-réduites en format long générés par 
#la fonction précédente. 
#L'année de départ sera la plus grande des premières années pour lesquelles on a des données.

#fun_data3_bis permet de prendre aussi en argument data_indiv (utile pour travailler sur les espèces)

fun_data3_bis<-function(data_list, #liste des tableaux longs contenant les données d'une covariable 
                    #pour toutes les années consécutives entre la première date de mesure 
                    #(coïncide avec une année d'inventaire) et la dernière (2016)
                    data_indiv=data_indiv #tableau des individus contenant les colonnes mort, recrutement, censure
){
  
  P <-
    length(data_list) + 1
  start <-
    max(unlist(lapply(data_list, function(x) {
      return(min(x$annee))
    })))
  #vérifier que la dernière année est la dernière année d'inventaire pour tout le monde
  if (any(unlist(lapply(data_list,
                        function(x) {
                          return(max(x$annee) != max(Dates_inv))
                        }))))
    stop("For all data the last year of measurements should be the last inventory : ",
         max(Dates_inv))
  #sélectionner les individus pour lesquels on a toutes les mesures
  ID_indiv <-
    unique(data_list[[1]]$ID)
  for (i in 1:length(data_list)) {
    ID_indiv <- intersect(ID_indiv, unique(data_list[[i]]$ID))
  }
  I <-
    length(ID_indiv)
  for (i in 1:length(data_list)) {
    data_list[[i]] <- data_list[[i]] %>% filter(ID %in% ID_indiv)
  }
  L <-
    max(Dates_inv) - start + 1
  
  Dates_start <-
    Dates_inv[which(Dates_inv >= start)]
  data_indiv_start <-
    data_indiv %>%
    filter(mort >
             start) %>%
    filter(ID %in% ID_indiv) %>%
    mutate(recrutement =
             sapply(recrutement, function(x)
               max(x, start))) %>%
    mutate(index_recr =
             sapply(recrutement, function(x)
               which(Dates_start == x))) %>%
    mutate(index_mort =
             sapply(mort, function(x)
               which(Dates_start == x)))
  
  #on modifie data_list pour ne prendre que les années supérieures à start 
  #et les individus dans data_indiv_start
  data_list <-
    lapply(data_list,
           function(x) {
             x %>%
               filter(x$annee >=
                        start) %>%
               spread(key =
                        annee, value = data) %>%
               filter(ID %in% data_indiv_start$ID) %>%  #
               arrange(ID)
           })
  I<-nrow(data_indiv_start) #
  X <-
    array(dim = c(I, L, P))
  for (i in 1:I) {
    for (j in 1:L) {
      X[i, j, 1] <- 1
      for (c in 2:P) {
        X[i, j, c] <-
          as.numeric(data_list[[c - 1]][i, j + 1]) #c-1 car la première covariable est l'intercept, j+1 car la première colonne est ID
      }
    }
  }
  
  
  return(
    list(
      D = length(Dates_start),
      Dates =
        Dates_start - (start - 1),
      I =
        nrow(data_indiv_start),
      P =
        dim(X)[3],
      L =
        dim(X)[2],
      X =
        X,
      Recr = data_indiv_start$index_recr,
      Mort = data_indiv_start$index_mort,
      Cens = data_indiv_start$censure
    ))
}

fun_data3<-function(data_list){
  return(
    fun_data3_bis(data_list, data_indiv)
  )}

## fun_data4 --------------------

#fun_data4 est l'exacte analogue de fun_data3 mais pour prendre aussi en compte 
#les données de censure de dépassement (arbres censurés lorsqu'ils dépassent 10cm de diamètre)
#(on remplace juste mort par mort2, censure par censure2 etc)

fun_data4<-function(data_list #liste des tableaux longs contenant les données d'une covariable 
                    #pour toutes les années consécutives entre la première date de mesure 
                    #(coïncide avec une année d'inventaire) et la dernière (2016)
){
  
  P <-
    length(data_list) + 1
  start <-
    max(unlist(lapply(data_list, function(x) {
      return(min(x$annee))
    })))
  #vérifier que la dernière année est la dernière année d'inventaire pour tout le monde
  if (any(unlist(lapply(data_list,
                        function(x) {
                          return(max(x$annee) != max(Dates_inv))
                        }))))
    stop("For all data the last year of measurements should be the last inventory : ",
         max(Dates_inv))
  #sélectionner les individus pour lesquels on a toutes les mesures
  ID_indiv <-
    unique(data_list[[1]]$ID)
  for (i in 1:length(data_list)) {
    ID_indiv <- intersect(ID_indiv, unique(data_list[[i]]$ID))
  }
  I <-
    length(ID_indiv)
  for (i in 1:length(data_list)) {
    data_list[[i]] <- data_list[[i]] %>% filter(ID %in% ID_indiv)
  }
  L <-
    max(Dates_inv) - start + 1
  
  Dates_start <-
    Dates_inv[which(Dates_inv >= start)]
  data_indiv_start <-
    data_indiv %>%
    filter(mort2 >
             start) %>%
    filter(ID %in% ID_indiv) %>%
    mutate(recrutement =
             sapply(recrutement, function(x)
               max(x, start))) %>%
    mutate(index_recr =
             sapply(recrutement, function(x)
               which(Dates_start == x))) %>%
    mutate(index_mort =
             sapply(mort2, function(x)
               which(Dates_start == x)))
  
  #on modifie data_list pour ne prendre que les années supérieures à start
  data_list <-
    lapply(data_list,
           function(x) {
             x %>%
               filter(x$annee >=
                        start) %>%
               spread(key =
                        annee, value = data) %>%
               arrange(ID)
           })
  X <-
    array(dim = c(I, L, P))
  for (i in 1:I) {
    for (j in 1:L) {
      X[i, j, 1] <- 1
      for (c in 2:P) {
        X[i, j, c] <-
          as.numeric(data_list[[c - 1]][i, j + 1]) #c-1 car la première covariable est l'intercept, j+1 car la première colonne est ID
      }
    }
  }
  
  
  return(
    list(
      D = length(Dates_start),
      Dates =
        Dates_start - (start - 1),
      I =
        nrow(data_indiv_start),
      P =
        dim(X)[3],
      L =
        dim(X)[2],
      X =
        X,
      Recr = data_indiv_start$index_recr,
      Mort = data_indiv_start$index_mort,
      Cens = data_indiv_start$censure2
    ))
}

## fun_data_diam ------------------------

# fun_data_diam est à utiliser avec script_gen_diam
# Il n'y a pas d'intercept (puisqu'il est compris dans le C0 du modèle)

fun_data_diam_bis<-function(data_list, #liste des tableaux longs contenant les données d'une covariable 
                        #pour toutes les années consécutives entre la première date de mesure 
                        #(coïncide avec une année d'inventaire) et la dernière (2016)
                        data_diam_inv = data_filled_diam_inv, #tableau long (ID, annee, data) 
                        #des inverses des diamètres, années consécutives
                        data_indiv=data_indiv #tableau des individus contenant les colonnes mort, recrutement, censure
){
  
  P <-
    length(data_list)
  start <-
    max(unlist(lapply(data_list, function(x) {
      return(min(x$annee))
    })))
  #vérifier que la dernière année est la dernière année d'inventaire pour tout le monde
  if (any(unlist(lapply(data_list,
                        function(x) {
                          return(max(x$annee) != max(Dates_inv))
                        }))))
    stop("For all data the last year of measurements should be the last inventory : ",
         max(Dates_inv))
  #sélectionner les individus pour lesquels on a toutes les mesures
  ID_indiv <-
    unique(data_list[[1]]$ID)
  for (i in 1:length(data_list)) {
    ID_indiv <- intersect(ID_indiv, unique(data_list[[i]]$ID))
  }
  I <-
    length(ID_indiv)
  for (i in 1:length(data_list)) {
    data_list[[i]] <- data_list[[i]] %>% filter(ID %in% ID_indiv)
  }
  L <-
    max(Dates_inv) - start + 1
  
  Dates_start <-
    Dates_inv[which(Dates_inv >= start)]
  data_indiv_start <-
    data_indiv %>%
    filter(mort >
             start) %>%
    filter(ID %in% ID_indiv) %>%
    mutate(recrutement =
             sapply(recrutement, function(x)
               max(x, start))) %>%
    mutate(index_recr =
             sapply(recrutement, function(x)
               which(Dates_start == x))) %>%
    mutate(index_mort =
             sapply(mort, function(x)
               which(Dates_start == x)))
  
  #on modifie data_list pour ne prendre que les années supérieures à start 
  #et les individus dans data_indiv_start
  data_list <-
    lapply(data_list,
           function(x) {
             x %>%
               filter(x$annee >=
                        start) %>%
               spread(key =
                        annee, value = data) %>%
               filter(ID %in% data_indiv_start$ID) %>%  #
               arrange(ID)
           })
  I<-nrow(data_indiv_start) #
  X <-
    array(dim = c(I, L, P))
  for (i in 1:I) {
    for (j in 1:L) {
      for (c in 1:P) {
        X[i, j, c] <-
          as.numeric(data_list[[c]][i, j + 1]) #j+1 car la première colonne est ID
      }
    }
  }
  
  Diam_inv<-data_diam_inv %>% 
    filter(ID %in% data_indiv_start$ID) %>% 
    filter(annee>=start) %>% 
    spread(key=annee, value=data) %>%
    arrange(ID) %>% 
    dplyr::select(-ID) %>% 
    as.matrix()
  
  return(
    list(
      D = length(Dates_start),
      Dates =
        Dates_start - (start - 1),
      I =
        nrow(data_indiv_start),
      P =
        dim(X)[3],
      L =
        dim(X)[2],
      X =
        X,
      Diam_inv = Diam_inv,
      Recr = data_indiv_start$index_recr,
      Mort = data_indiv_start$index_mort,
      Cens = data_indiv_start$censure
    ))
}

fun_data_diam<-function(data_list){
  return(
    fun_data_diam_bis(data_list, data_filled_diam_inv, data_indiv)
  )}


# Pour prendre en compte la censure par accident :
fun_data_diam_acc<-function(data_list){
  data_indiv_acc<-data_indiv %>% 
    mutate(mort=mort_acc, censure=censure_acc, index_mort=index_mort_acc)
  return(
    fun_data_diam_bis(data_list, 
                      (data_filled_diam_inv %>% filter(ID %in% data_indiv_acc$ID)), 
                      data_indiv_acc)
  )}

# Pour pouvoir filtrer par espèce :
fun_data_diam_acc_bis<-function(data_list, 
                                data_filled_diam_inv=data_filled_diam_inv,
                                data_indiv=data_indiv){
  data_indiv_acc<-data_indiv %>% 
    mutate(mort=mort_acc, censure=censure_acc, index_mort=index_mort_acc)
  return(
    fun_data_diam_bis(data_list, 
                      (data_filled_diam_inv %>% filter(ID %in% data_indiv_acc$ID)), 
                      data_indiv_acc)
  )
}

# Pour prendre en compte la censure par accident et la censure 2 (lorsque le DBH dépasse 10 cm)
fun_data_diam_acc_cens2_bis<-function(data_list, 
                                  data_filled_diam_inv=data_filled_diam_inv,
                                  data_indiv=data_indiv){
  data_indiv_temp<-data_indiv %>% 
    mutate(mort=pmin(mort_acc, mort2),
           censure=pmax(censure_acc, censure2),
           index_mort=pmin(index_mort_acc, index_mort2))
  return(
    fun_data_diam_bis(data_list, 
                      (data_filled_diam_inv %>% filter(ID %in% data_indiv_temp$ID)), 
                      data_indiv_temp)
  )  
}

fun_data_diam_acc_cens2<-function(data_list){
  fun_data_diam_acc_cens2_bis(data_list, data_filled_diam_inv, data_indiv)
}
