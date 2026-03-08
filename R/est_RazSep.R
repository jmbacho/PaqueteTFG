#' @title Estimador de razón separado en Muestreo Aleatorio Simple Estratificado (MASE)
#' @description Calcula las estimaciones de razón, media y total, así como sus varianzas estimadas utilizando el estimador de razón separado. Esta función es específica para diseños MASE.
#' @param y Vector muestral de la variable de interés. Los datos deben aparecer ordenados por estratos.
#' @param x Vector muestral de la variable auxiliar. Debe tener la misma longitud que \code{y} y estar ordenado por estratos.
#' @param nh Vector que indica el tamaño muestral de cada estrato.
#' @param Nh Vector que indica el tamaño poblacional de cada estrato.
#' @param tXh Vector que indica el total poblacional conocido de la variable auxiliar X en cada estrato.
#' @details
#' La función implementa el estimador de razón separado, realizando la estimación de la razón en cada estrato para posteriormente formar el promedio de estas estimaciones separadas como una sola estimación de la razón poblacional.
#' Como se ha indicado anteriormente, se asume que los vectores muestrales están ordenados por estratos según la secuencia el vector \code{nh}.
#' @return Una lista con los siguientes componentes:
#' \itemize{
#'    \item \code{Estimador de razón separado en cada estrato}
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
#' tXh <- c(11,24)
#' #datos muestrales
#' nh <- c(2,2)
#' #vectores muestrales ordenados por estrato
#' y <- c(1,2,4,5)
#' x <- c(2,4,5,7)
#' est_RazSep(y,x,nh,Nh,tXh)

est_RazSep <- function(y,x,nh,Nh,tXh)
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
  if(length(Nh)!=L || length(tXh)!=L) stop("Error, la dimensión de los parámetros 'Nh' y 'tXh' deben ser la misma.")
  
  #vector que indica a qué estrato pertenece cada observación muestral
  ind <- rep(1:L,nh)
  
  #extraemos los valores muestrales por estrato
  xh <- split(x,ind)
  yh <- split(y,ind)
  
  #calculamos los estimadores de razón en cada estrato
  estRh <- mapply(function(y,x) mean(y)/mean(x),yh,xh)
  
  #calculamos los pesos en cada estrato
  Wh <- Nh/sum(Nh)
  
  #calculamos la estimación de la media mediante el estimador de razón separado
  estRsMedia <- sum(Wh*estRh*tXh/Nh)
  estRsTotal <- estRsMedia*sum(Nh)
  
  #calculamos las fracciones de muestreo en cada estrato
  fh <- nh/Nh
  
  #calculamos las estimaciones de la varianza correspondientes
  var_yh <- sapply(yh,var)
  var_xh <- sapply(xh,var)
  cov_yxh <- mapply(cov,yh,xh)
  varEstMedia <- sum(Wh^2 * (1-fh)/nh * (var_yh + estRh^2*var_xh - 2*estRh*cov_yxh))
  varEstTotal <- varEstMedia*sum(Nh)^2
  
  #devolvemos los valores
  return(list("Estimador de razón separado en cada estrato"=estRh,
              "Estimación de la media"=estRsMedia,
              "Estimación de la varianza media"=varEstMedia,
              "Estimación del total"=estRsTotal,
              "Estimación de la varianza total"=varEstTotal))
}