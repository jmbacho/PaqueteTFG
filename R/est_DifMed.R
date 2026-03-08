#' @title Estimador de diferencia de medias en Muestreo Aleatorio Simple (MAS)
#' @description Calcula las estimaciones de la media y el total poblacional, así como sus varianzas estimadas, utilizando el estimador de diferencia. Este es un caso particular del estimador de regresión donde se asume una pendiente fija b=1.
#' @param y Vector con los valores de la muestra observada de la variable de interés.
#' @param x Vector con los valores de la muestra de la variable auxiliar. Debe tener la misma longitud que \code{y}.
#' @param tX Valor que indica el total poblacional conocido de la variable auxiliar X.
#' @param N Valor que indica el tamaño total de la población.
#' @details
#' Es eficiente cuando la variable auxiliar X y la variable de interés Y tienen una correlación muy alta y una relación de pendiente igual a 1
#' @return Una lista con los siguientes componentes:
#' \itemize{
#'    \item \code{Estimación de la media}
#'    \item \code{Estimación de la varianza media}
#'    \item \code{Estimación del total}
#'    \item \code{Estimación de la varianza total}
#' }
#' @author Javier Miralles Corbacho
#' @export
#' @examples
#' #datos poblacionales 
#' N <- 100
#' tX <- 500
#' #datos muestrales (x e y muy parecidos, luego b approx. 1)
#' y <- c(51,49,52,50,48)
#' x <- c(50,48,51,50,49)
#' est_DifMed(y,x,tX,N)

est_DifMed <- function(y,x,tX,N)
{
  #validaciones de los parámetros sobre:
  #1.- los vectores muestrales 'y' y 'x': han de ser vectores, han de tener la misma longitud y no deben contener NA valores
  if((is.vector(y) && is.vector(x))==FALSE) stop("Error, los parámetros 'y' y 'x' han de ser vectores.")
  if(length(y)!=length(x)) stop("Error, la longitud de los vectores muestrales ha de ser la misma.")
  if(any(is.na(y)) || any(is.na(x))) stop("Error, los vectores no deben contener valores NA.")
  
  #2.- el tamaño poblacional 'N': que sea mayor que el tamaño muestral y que sea un valor mayor o igual que 2
  if((N<1)==TRUE) stop("Error, el tamaño poblacional ha de ser un número positivo mayor que 2.")
  if(length(y)>N) stop("Error, el tamaño muestral ha de ser inferior al tamaño poblacional 'N'.")
  
  #3.- el total poblacional 'tX': ha de ser un número mayor que 1
  if(is.numeric(tX)==FALSE || length(tX)!=1 || tX<=1) stop("Error, el parámetro 'tX' ha de ser un número mayor que 1.")
  
  #calculamos el tamaño muestral
  n <- length(y)
  
  #calculamos la fracción de muestreo
  f <- n/N
  
  #calculamos la estimación de la media y del total
  estMedia <- mean(y) + (tX/N-mean(x))
  estTotal <- estMedia*N
  
  #calculamos las estimaciones de la varianza correspondientes
  varEstMedia <- (1-f)/n * (var(y)+var(x)-2*cov(x,y))
  varEstTotal <- varEstMedia*N^2
  
  #devolvemos los valores
  return(list("Estimación de la media"=estMedia,
              "Estimación de la varianza media"=varEstMedia,
              "Estimación del total"=estTotal,
              "Estimación de la varianza total"=varEstTotal
  ))
}