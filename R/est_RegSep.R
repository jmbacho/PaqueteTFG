#' @title Estimador de regresión separado en Muestreo Aleatorio Simple Estratificado (MASE)
#' @description Calcula las estimaciones de la pendiente de regresión, media y total, así como sus varianzas estimadas utilizando el estimador de regresión separado. Esta función es específica para diseños MASE.
#' @param y Vector muestral de la variable de interés. Los datos deben aparecer ordenados por estratos.
#' @param x Vector muestral de la variable auxiliar. Debe tener la misma longitud que \code{y} y estar ordenado por estratos.
#' @param nh Vector que indica el tamaño muestral de cada estrato.
#' @param Nh Vector que indica el tamaño poblacional de cada estrato.
#' @param tXh Vector que indica el total poblacional conocido de la variable auxiliar X en cada estrato.
#' @details
#' La función implementa el estimador de regresión separado, realizando la estimación de la pendiente de regresión en cada estrato para posteriormente agregar las estimaciones de media y total ponderadas por el tamaño de los estratos.
#' @return Una lista con los siguientes componentes:
#' \itemize{
#'    \item \code{Pendiente de regresión separada en cada estrato}
#'    \item \code{Estimación de la media}
#'    \item \code{Estimación de la varianza media}
#'    \item \code{Estimación del total}
#'    \item \code{Estimación de la varianza total}
#' }
#' @author Javier Miralles Corbacho
#' @export
#' @examples
#' #datos poblacionales (2 estratos)
#' Nh <- c(400, 350)      #tamaño poblacional por estrato
#' tXh <- c(2000, 1800)   #total poblacional de X por estrato
#' #datos muestrales (4 observaciones del primer estrato y 3 del segundo)
#' nh <- c(4, 3)
#' #vectores muestrales ordenados por estrato
#' y <- c(12,14,15,13,20,22,19)
#' x <- c(4.5,5.0,5.5,4.8,8.0,8.5,7.9)
#' est_RegSep(y,x,nh,Nh,tXh)

est_RegSep <- function(y,x,nh,Nh,tXh)
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
  if(length(y)!=sum(nh)) stop("Error, las dimensiones de los vectores muestrales 'x' e 'y' no coincide con la dimensión indicada por el vector de tamaños 'nh'.")
  
  #3.- el vector de tamaño poblacional de los estratos 'Nh': ha de ser mayor que nh valor a valor
  if(any(Nh<nh)) stop("Error, el tamaño poblacional no puede ser menor que el tamaño muestral.")
  if(length(Nh)!=L || length(tXh)!=L) stop("Error, la dimensión de los parámetros 'Nh' y 'tXh' deben ser la misma.")
  
  #vector que indica a qué estrato pertenece cada observación muestral
  ind <- rep(1:L,nh)
  
  #extraemos los valores muestrales por estrato
  xh <- split(x,ind)
  yh <- split(y,ind)
  
  #calculamos los estimadores de razón en cada estrato, comprobando que se pueda calcular la estimación de la pendiente de la recta de regresión en cada estrato
  var_xh <- sapply(xh,var)
  #determinamos el índice del estrato en el que se ha encontrado (en caso de que haya) un vector xh cuya varianza sea cero (es contante)
  indErr <- which(var_xh==0)[1]
  if(any(var_xh==0)) stop(paste0("Error, no se puede calcular la estimación de la pendiente de la recta de regresión en el estrato número ", indErr, " debido a que el vector de parámetros en dicho estrato es constante (tiene varianza cero)."))
  estbh <- mapply(function(y,x) cov(y,x)/var(x),yh,xh)
  
  #calculamos los pesos en cada estrato
  Wh <- Nh/sum(Nh)
  
  #calculamos la media de los vectores x e y en cada estrato
  m_xh <- sapply(xh,mean)
  m_yh <- sapply(yh,mean)
  #calculamos la media poblacional de X en cada estrato
  med_Xh <- tXh/Nh
  #calculamos la estimación de regresión en cada estrato
  estMedia_h <- m_yh + estbh*(med_Xh-m_xh)
  #calculamos la estimación de la media realizando la suma ponderada
  estRegsMedia <- sum(Wh*estMedia_h)
  
  #calculamos la estimación del total
  estRegsTotal <- estRegsMedia*sum(Nh)
  
  #calculamos las fracciones de muestreo en cada estrato
  fh <- nh/Nh
  
  #calculamos las estimaciones de la varianza correspondientes
  var_yh <- sapply(yh,var)
  var_xh <- sapply(xh,var)
  cov_yxh <- mapply(cov,yh,xh)
  varEstMedia <- sum(Wh^2 * (1-fh)/nh * (var_yh + estbh^2*var_xh - 2*estbh*cov_yxh))
  varEstTotal <- varEstMedia*sum(Nh)^2
  
  #devolvemos los valores
  return(list("Estimador de pendiente de regresión separada en cada estrato"=estbh,
              "Estimación de la media"=estRegsMedia,
              "Estimación de la varianza media"=varEstMedia,
              "Estimación del total"=estRegsTotal,
              "Estimación de la varianza total"=varEstTotal))
}
