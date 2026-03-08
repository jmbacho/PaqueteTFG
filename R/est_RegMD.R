#' @title Estimador de Regresión en Muestreo Doble
#' @description Calcula las estimaciones de la media y el total poblacional, así como sus varianzas estimadas, utilizando un estimador de regresión lineal en un diseño de muestreo doble.
#' @param xp Vector muestral de la variable auxiliar (obtenido en la primera fase).
#' @param y Vector muestral de la variable de interés (obtenido en la segunda fase).
#' @param x Vector muestral de la variable auxiliar (obtenido en la segunda fase). Debe tener la misma longitud que \code{y}.
#' @param N Valor que indica el tamaño total de la población.
#' @details
#' En el muestreo doble para regresión, se utiliza una primera muestra para estimar la media poblacional de la variable auxiliar, y una segunda muestra más pequeña para estimar la pendiente de regresión y ajustar la media de la variable de interés.
#' @return Una lista con los siguientes componentes:
#' \itemize{
#'    \item \code{Estimación de la pendiente (b)}
#'    \item \code{Estimación de la media de X (fase 1)}
#'    \item \code{Estimación de la media}
#'    \item \code{Estimación de la varianza media}
#'    \item \code{Estimación del total}
#'    \item \code{Estimación de la varianza total}
#'    \item \code{Coeficiente de correlación de Pearson entre x e y}
#' }
#' @author Javier Miralles Corbacho
#' @export
#' @examples
#' #datos poblacionales
#' N <- 1000
#' #fase1 (muestra grande de la variable auxiliar x, n'=100)
#' xp <- rnorm(100,50,5)
#' #fase2 (submuestra de tamaño 20)
#' x <- xp[1:20] 
#' y <- 2*x + rnorm(20,0,1)
#' est_RegMD(xp,y,x,N)

est_RegMD <- function(xp,y,x,N)
{
  #validaciones de los parámetros sobre:
  #1.- los vectores muestrales 'y', 'x' y 'xp': han de ser vectores y no contener NA valores
  if((is.vector(y) && is.vector(x) && is.vector(xp))==FALSE) stop("Error, los parámetros 'y', 'x' y 'xp' han de ser vectores.")
  if(any(is.na(y)) || any(is.na(x)) || any(is.na(xp))) stop("Error, los vectores no deben contener NA valores.")
  
  #2.- vectores 'x' e 'y': deben tener la misma longitud
  if(length(y) != length(x)) stop("Error, los vectores de segunda fase 'y' y 'x' deben tener la misma longitud.")
  
  #3.- tamaño muestral de fase 1 'xp' y tamaño muestral de fase 2 'n': el tamaño muestral de la fase dos ha de ser inferior al de la fase 1 y ambos inferiores al tamaño poblacional
  if(length(y) >= length(xp)) stop("Error, el tamaño muestral de segunda fase 'n' debe ser inferior al de primera fase 'np'.")
  if(length(xp) >= N) stop("Error, el tamaño de la muestra de primera fase 'np' ha de ser menor que el tamaño poblacional N.")
  
  #fase1
  #calculamos la media de la variable auxiliar
  xpmed <- mean(xp)
  #determinamos el tamaño muestral de la variable auxiliar
  np <- length(xp)
  
  #fase2
  #determinamos el tamaño muestral de la variable de interés
  n <- length(y)
  #calculamos el estimador de la pendiente de la recta de regresión
  if(var(x)==0) stop("Error, no se puede calular la estimación de la pendiente de la recta de regresión debido a que el vector de parámetros 'x' es constante (tiene varianza cero).")
  estb <- cov(x,y)/var(x)
  #calculamos el coeficiente de correlación de Pearson entre X e Y
  rho <- cor(x,y)
  
  #calculamos ahora ambas estimaciones
  estMedia <- mean(y) + estb*(xpmed-mean(x))
  estTotal <- estMedia*N
  
  #calculamos las estimaciones de la varianza correspondientes
  varEstMedia <- (N-np)/(N*np)*var(y) + (np-n)/(np*n)*(1-rho^2)*var(y)
  varEstTotal <- varEstMedia*N^2
  
  #devolvemos los valores
  return(list("Estimación de la pendiente"=estb,
              "Estimación de la media de X (fase 1)"=xpmed,
              "Estimación de la media"=estMedia,
              "Estimación de la varianza media"=varEstMedia,
              "Estimación del total"=estTotal,
              "Estimación de la varianza total"=varEstTotal,
              "Coeficiente de correlación de Pearson entre x e y"=rho
  ))
}