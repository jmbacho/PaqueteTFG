#' @title Estimador de regresión en Muestreo Aleatorio Simple (MAS)
#' @description Calcula las estimaciones de la media y el total poblacional, así como sus varianzas estimadas, utilizando el estimador de regresión lineal.
#' @param y Vector con los valores de la muestra observada de la variable de interés.
#' @param x Vector con los valores de la muestra de la variable auxiliar. Debe tener la misma longitud que \code{y}.
#' @param tX Valor que indica el total poblacional conocido de la variable auxiliar X.
#' @param N Valor que indica el tamaño total de la población.
#' @details
#' La función calcula primero la pendiente óptima de regresión basada en la covarianza muestral para posteriormente aplicar la fórmula del estimador de regresión.
#' Se trata de un estimador es más eficiente que el de razón cuando la recta de regresión entre la variable de interés y la auxiliar no pasa por el origen.
#' @return Una lista con los siguientes componentes:
#' \itemize{
#'    \item \code{Pendiente estimada}
#'    \item \code{Estimación de la media}
#'    \item \code{Estimación de la varianza media}
#'    \item \code{Estimación del total}
#'    \item \code{Estimación de la varianza total}
#'    \item \code{Coeficiente de correlación de Pearson entre x e y}
#' }
#' @author Javier Miralles Corbacho
#' @export
#' @examples
#' #datos poblacionales (se quiere estimar ganancias basándose (Y) en las del año anterior (X))
#' N <- 123      #total de franquicias
#' tX <- 128200  #total ventas año anterior
#' #datos muestrales (13 franquicias)
#' x <- c(550,720,1500,1020,620,980,928,1200,1350,1750,670,728,1530)
#' y <- c(610,780,1600,1030,600,1050,977,1440,1570,2210,980,865,1710)
#' est_RegMas(y,x,tX,N)

est_RegMas <- function(y,x,tX,N)
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
  
  #calculamos el coeficiente de correlación de Pearson entre X e Y
  rho <- cor(x,y)
  
  #calculamos el estimador de la pendiente de la recta de regresión
  if(var(x)==0) stop("Error, no se puede calcular la estimación de la pendiente de la recta de regresión debido a que el vector de parámetros 'x' es constante (tiene varianza cero).")
  estb <- cov(x,y)/var(x)
  
  #calculamos la estimación de la media y del total mediante el estimador de regresión
  estMedia <- mean(y) + estb*(tX/N-mean(x))
  estTotal <- estMedia*N
  
  #calculamos las estimaciones de la varianza correspondientes
  varEstMedia <- (1-f)/n * (var(y)+estb^2*var(x)-2*estb*cov(x,y))
  varEstTotal <- varEstMedia*N^2
  
  #devolvemos los valores
  return(list("Estimador de la pendiente"=estb,
              "Estimación de la media"=estMedia,
              "Estimación de la varianza media"=varEstMedia,
              "Estimación del total"=estTotal,
              "Estimación de la varianza total"=varEstTotal,
              "Coeficiente de correlación de Pearson entre x e y"=rho
  ))
}