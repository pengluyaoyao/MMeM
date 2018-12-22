#' parses formulas to creates model matrices
#' @param fml a two-sided linear formula object describing both the fixed-effects and random-effects  parts  of  the  model,
#'  with  the  response  on  the  left  of  a ~ operator. For univariate response, put variable name directly; for multivariate responses
#'  combine variables using concatenate operator, for example, for bivariate responses, c(var1, var2).  The predictor terms are separated  by + operators,  on  the  right.   Random-effects  terms  are
#' distinguished by vertical bars '|' separating expressions for design matrices from grouping factors.
#' @param data data frame containing the variables named in formula.
#' @param factor_X (logical) indicating whether predictor is a factor or continuous. By default is TRUE
#' @importFrom lme4 findbars
#' @importFrom stringr str_extract_all
#' @export
MMeM_terms <- function(fml, data, factor_X){

  Fml = stats::formula(fml)
  any_RE <- length( lme4::findbars(Fml))

  if(any_RE == 0){
    stop('No random effects in the model')
  } else{
    df = stats::get_all_vars(Fml, data)

    terms = attr(stats::terms.formula(Fml), 'variables')
    DVs = all.names(terms[2])
    if(length(DVs) == 3){
      message(paste('Bivariate response:', DVs[2], 'and', DVs[3]))
      Y = as.matrix(df[,match(DVs[2:length(DVs)], colnames(df))])
      DV = c(DVs[2], DVs[3])
    }else if(length(DVs) == 1){
      message(paste('Univariate response:', DVs[1:length(DVs)]))
      Y = as.matrix(df[,match(DVs[1:length(DVs)], colnames(df))])
      DV = DVs[1]
    }else{
      stop('Dependent variables should be univariate or bivariate.')
    }

    re_term = stringr::str_extract_all(format(terms[length(terms)]), "(?<=\\|).+?(?=\\))")[[1]]
    re_strip = gsub(" ", "", re_term, fixed = TRUE)
    re_data = df[,match(re_strip, colnames(df))]
    Z = stats::model.matrix(~ -1 + factor(re_data))

    IV = all.names(terms[-c(1,2,length(terms))])
    IV_data = df[,match(IV, colnames(df))]
    if(factor_X == TRUE){
      X = stats::model.matrix(~ factor(IV_data))
    }else{
      X = stats::model.matrix(~ IV_data)
    }


    N = nrow(X)
    I = diag(1, N)
    q = ncol(as.matrix(Y))
  }

  return(list(X = X, Z = Z, Y = Y, N =N, I =I, q=q, DV = DV))
}


#d = MMeM_terms(c(age,male) ~ male+(1|id), data=alcohol1)



