#' @title Estimador de Muestreo Doble para Estratificación (MASE)
#' @description Calcula las estimaciones de la media y el total poblacional, así como sus varianzas estimadas, utilizando un diseño de muestreo doble (bifásico) para estratificación.
#' @param y Vector muestral de la variable de interés (obtenido en la segunda fase). Los datos deben aparecer ordenados por estratos.
#' @param nh Vector que indica el tamaño muestral de cada estrato en la segunda fase.
#' @param nhp Vector que indica el tamaño muestral de cada estrato en la primera fase (n').
#' @param N Valor numérico que indica el tamaño total de la población.
#' @details
#' En el muestreo doble para estratificación, se utiliza una primera muestra grande para estimar los pesos de los estratos, y una segunda muestra más pequeña para estimar la media de la variable de interés en cada estrato.
#' @return Una lista con los siguientes componentes:
#' \itemize{
#'    \item \code{Pesos estimados de los estratos (fase 1)}
#'    \item \code{Estimación de la media}
#'    \item \code{Estimación de la varianza media}
#'    \item \code{Estimación del total}
#'    \item \code{Estimación de la varianza total}
#' }
#' @author Javier Miralles Corbacho
#' @export
#' @examples
#' #datos poblacionales (2 estratos)
#' N <- 10000
#' #fase1 (muestra de tamaño 100, 60 del primer estrato y 40 del segundo)
#' nhp <- c(60,40)
#' #fase2 (submuestra de tamaño 20, 10 del primer estrato y 10 del segundo)
#' y <- c(10,12,10,12,10,12,10,12,10,12,20,24,20,24,20,24,20,24,20,24)
#' nh <- c(10,10)
#' est_EstrMD(y,nh,nhp,N)

est_EstrMD <- function(y,nh,nhp,N)
{
  #validaciones de los parámetros sobre:
  #1.- el vector de muestra 'y': ha de ser vector y no debe contener NA valores
  if(is.vector(y)==FALSE) stop("Error, el parámetro 'y' ha de ser un vector.")
  if(any(is.na(y))==TRUE) stop("Error, el vector de muestra 'y' contiene NA valores.")
  
  #2.- los vectores de tamaños muestrales 'nh' y 'nhp': deben tener la misma longitud y el número de estratos ha de ser mayor o igual que 2
  if(length(nh) != length(nhp)) stop("Error, los vectores 'nh' y 'nhp' deben tener la misma longitud (mismo número de estratos).")
  L <- length(nh) # Número de estratos
  if(L<=1) stop("Error, el número de estratos ha de ser un número positivo natural mayor o igual que 2.")
  
  #3.- tamaño muestral
  if(sum(nh) != length(y)) stop("Error, la dimensión del vector 'y' no coincide con la suma de los tamaños muestrales 'nh'.")
  if(any(nh>nhp)) stop("Error, el tamaño muestral de segunda fase 'nh' no puede ser superior al de primera fase 'nhp' en ningún estrato.")
  if(sum(nhp)>=N) stop("Error, el tamaño de la muestra de primera fase ha de ser menor que el tamaño poblacional N.")
  
  #fase1
  #calculamos el tamaño muestral
  np <- sum(nhp)
  
  #calculamos los pesos estimados en cada estrato
  Whp <- nhp/np
  
  #fase2
  #vector que indica a qué estrato pertenece cada observación muestral
  ind <- rep(1:L, nh)
  #extraemos los valores muestrales por estrato
  yh <- split(y,ind)
  
  #calculamos la media muestral en cada estrato
  yhmed <- sapply(yh,mean)

  #calculamos la estimación de la media y del total
  estMedia <- sum(Whp*yhmed)
  estTotal <- estMedia*N
  
  #calculamos las estimaciones de la varianza correspondientes
  varEstMedia <- (N-np)/(N*(np-1)) * sum(Whp*(yhmed-estMedia)^2) + (N-1)/N * sum(((nhp-1)/(np-1) - (nh-1)/(N-1))*Whp*sapply(yh,var)/nh)
  varEstTotal <- varEstMedia*N^2
  
  #devolvemos los valores
  return(list("Pesos estimados de los estratos (fase 1)"=Whp,
              "Estimación de la media"=estMedia,
              "Estimación de la varianza media"=varEstMedia,
              "Estimación del total"=estTotal,
              "Estimación de la varianza total"=varEstTotal
  ))
}