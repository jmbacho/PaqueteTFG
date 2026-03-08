#' @title Estimador de razón combinado en Muestreo Aleatorio Simple Estratificado (MASE)
#' @description Calcula las estimaciones de razón, media y total, así como sus varianzas estimadas utilizando el estimador de razón combinado.
#' @param y Vector muestral de la variable de interés. Los datos deben aparecer ordenados por estratos.
#' @param x Vector muestral de la variable auxiliar. Debe tener la misma longitud que \code{y} y estar ordenado por estratos.
#' @param nh Vector que indica el tamaño muestral de cada estrato.
#' @param Nh Vector que indica el tamaño poblacional de cada estrato.
#' @param tX Número que indica el total poblacional conocido de la variable auxiliar X.
#' @details
#' La función calcula primero las medias ponderadas globales para obtener una razón combinada, para posteriormente estimar el total y la media poblacional.
#' @return Una lista con los siguientes componentes:
#' \itemize{
#'    \item \code{Estimador de razón combinado}
#'    \item \code{Estimación de la media}
#'    \item \code{Estimación de la varianza media}
#'    \item \code{Estimación del total}
#'    \item \code{Estimación de la varianza total}
#' }
#' @author Javier Miralles Corbacho
#' @export
#' @examples
#' #datos poblacionales (2 estratos)
#' Nh <- c(3,3)     
#' tX <- 35
#' #datos muestrales
#' nh <- c(2,2)
#' #vectores muestrales ordenados por estrato
#' y <- c(1,2,4,5)
#' x <- c(2,4,5,7)
#' est_RazComb(y,x,nh,Nh,tX)

est_RazComb <- function(y,x,nh,Nh,tX)
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
  
  #calculamos el estimador de razón combinado para la media
  estRc <- estY/estX
  
  #calculamos la estimación de la media y del total mediante el estimador de razón combinado
  estRcMedia <- estRc*tX/sum(Nh)
  estRcTotal <- estRc*tX
  
  #calculamos las fracciones de muestreo en cada estrato
  fh <- nh/Nh
  
  #calculamos las estimaciones de la varianza correspondientes
  var_yh <- sapply(yh,var)
  var_xh <- sapply(xh,var)
  cov_yxh <- mapply(cov,yh,xh)
  varEstMedia <- sum(Wh^2 * (1-fh)/nh * (var_yh + estRc^2*var_xh - 2*estRc*cov_yxh))
  varEstTotal <- varEstMedia*sum(Nh)^2
  
  #devolvemos los valores
  return(list("Estimador de razón combinado"=estRc,
              "Estimación de la media"=estRcMedia,
              "Estimación de la varianza media"=varEstMedia,
              "Estimación del total"=estRcTotal,
              "Estimación de la varianza total"=varEstTotal
  ))
}