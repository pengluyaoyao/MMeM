
MMeM_terms <- function(fml, data, factor_X){

  Fml = formula(fml)
  any_RE <- length( lme4::findbars(Fml))

  if(any_RE == 0){
    stop('No random effects in the model')
  } else{
    df = get_all_vars(Fml, data)

    terms = attr(terms.formula(Fml), 'variables')
    DVs = all.names(terms[2])
    if(length(DVs) == 3){
      message(paste('Bivariate response:', DVs[2], 'and', DVs[3]))
      Y = as.matrix(df[,match(DVs[2:length(DVs)], colnames(df))])
    }else if(length(DVs) == 1){
      message(paste('Univariate response:', DVs[1:length(DVs)]))
      Y = as.matrix(df[,match(DVs[1:length(DVs)], colnames(df))])
    }else{
      stop('Dependent variables should be univariate or bivariate.')
    }

    re_term = str_extract_all(format(terms[length(terms)]), "(?<=\\|).+?(?=\\))")[[1]]
    re_strip = gsub(" ", "", re_term, fixed = TRUE)
    re_data = df[,match(re_strip, colnames(df))]
    Z = model.matrix(~ -1 + factor(re_data))

    IV = all.names(terms[-c(1,2,length(terms))])
    IV_data = df[,match(IV, colnames(df))]
    if(factor_X == TRUE){
      X = model.matrix(~ factor(IV_data))
    }else{
      X = model.matrix(~ IV_data)
    }


    N = nrow(X)
    I = diag(1, N)
    q = ncol(as.matrix(Y))
  }

  return(list(X = X, Z = Z, Y = Y, N =N, I =I, q=q))
}


#d = MMeM_terms(c(age,male) ~ male+(1|id), data=alcohol1)



