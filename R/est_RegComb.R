#' @title Estimador de regresión combinado en Muestreo Aleatorio Simple Estratificado (MASE)
#' @description Calcula las estimaciones de la pendiente de regresión, media y total, así como sus varianzas estimadas utilizando el estimador de regresión combinado.
#' @param y Vector muestral de la variable de interés. Los datos deben aparecer ordenados por estratos.
#' @param x Vector muestral de la variable auxiliar. Debe tener la misma longitud que \code{y} y estar ordenado por estratos.
#' @param nh Vector que indica el tamaño muestral de cada estrato.
#' @param Nh Vector que indica el tamaño poblacional de cada estrato.
#' @param tX Número que indica el total poblacional conocido de la variable auxiliar X.
#' @details
#' La función calcula primero la pendiente de regresión combinada utilizando las varianzas y covarianzas ponderadas de todos los estratos, para posteriormente estimar el total y la media poblacional mediante la fórmula de regresión lineal.
#' @return Una lista con los siguientes componentes:
#' \itemize{
#'    \item \code{Pendiente de regresión combinada}
#'    \item \code{Estimación de la media}
#'    \item \code{Estimación de la varianza media}
#'    \item \code{Estimación del total}
#'    \item \code{Estimación de la varianza total}
#' }
#' @author Javier Miralles Corbacho
#' @export
#' @examples
#' #datos poblacionales (2 estratos)
#' Nh <- c(400, 350)  #tamaño poblacional por estrato
#' tX <- 3800         #total poblacional de X 
#' #datos muestrales (4 observaciones del primer estrato y 3 del segundo)
#' nh <- c(4, 3)
#' #vectores muestrales ordenados por estrato
#' y <- c(12,14,15,13,20,22,19)
#' x <- c(4.5,5.0,5.5,4.8,8.0,8.5,7.9)
#' est_RegComb(y,x,nh,Nh,tX)

est_RegComb <- function(y,x,nh,Nh,tX)
{
  #validaciones de los parámetros sobre:
  #1.- los vectores muestrales 'y' y 'x': han de ser vectores, han de tener la misma longitud y no deben contener NA valores
  if((is.vector(y) && is.vector(x))==FALSE) stop("Error, los parámetros 'y' y 'x' han de ser vectores.")
  if(length(y)!=length(x)) stop("Error, la longitud de los vectores muestrales ha de ser la misma.")
  if(any(is.na(y)) || any(is.na(x))) stop("Error, los vectores no deben contener valores NA.")
  
  #2.- el vector de tamaño muestral de los estratos 'nh': ha de tener longitud superior a 1 y han de sumar la longitud del vector muestral
  #definimos previamente el número de estratos
  L <- length(nh)
  #ahora procedemos con las validaciones
  if(L<=1) stop("Error, el número de estratos ha de ser un número positivo natural mayor o igual que 2.")
  if(length(y)!=sum(nh)) stop("Error, las dimensiones de los vectores muestrales 'x' e 'y' no conincide con la dimensión indicada por el vector de tamaños 'nh'.")
  
  #3.- el vector de tamaño poblacional de los estratos 'Nh': ha de ser mayor que nh valor a valor
  if(any(Nh<nh)) stop("Error, el tamaño poblacional no puede ser menor que el tamaño muestral.")
  
  #4.- el total poblacional 'tX': ha de ser un número mayor que 1
  if(is.numeric(tX)==FALSE || length(tX)!=1 || tX<=1) stop("Error, el parámetro 'tX' ha de ser un número mayor que 1.")
  
  #vector que indica a qué estrato pertenece cada observación muestral
  ind <- rep(1:L,nh)
  
  #extraemos los valores muestrales por estrato
  xh <- split(x,ind)
  yh <- split(y,ind)
  
  #calculamos los pesos en cada estrato
  Wh <- Nh/sum(Nh)
  
  #calculamos la media muestral ponderada de Y y de X
  estY <- sum(Wh*sapply(yh,mean))
  estX <- sum(Wh*sapply(xh,mean))
  
  #calculamos las fracciones de muestreo en cada estrato
  fh <- nh/Nh
  
  #calculamos el estimador de la pendiente de regresión combinado
  var_yh <- sapply(yh,var)
  var_xh <- sapply(xh,var)
  cov_yxh <- mapply(cov,yh,xh)
  estb <- sum(Wh^2*(1-fh)/nh*cov_yxh)/sum(Wh^2*(1-fh)/nh*var_xh)
    
  #calculamos la estimación de la media y del total mediante el estimador de regresión combinado
  estRegcMedia <- estY + estb*(tX/sum(Nh)-estX)
  estRegcTotal <- estRegcMedia*sum(Nh)
  
  #calculamos las estimaciones de la varianza correspondientes
  varEstMedia <- sum(Wh^2 * (1-fh)/nh * (var_yh + estb^2*var_xh - 2*estb*cov_yxh))
  varEstTotal <- varEstMedia*sum(Nh)^2
  
  #devolvemos los valores
  return(list("Estimador de regresión combinado"=estb,
              "Estimación de la media"=estRegcMedia,
              "Estimación de la varianza media"=varEstMedia,
              "Estimación del total"=estRegcTotal,
              "Estimación de la varianza total"=varEstTotal
  ))
}