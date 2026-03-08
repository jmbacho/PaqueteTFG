#' @title Estimador de Hájek para la media y el total
#' @description Calcula el estimador de Hajek para la media y el total poblacional.
#' @param y Vector con los valores de la muestra observada. No debe contener NAs.
#' @param N Valor que indica el tamaño de la población.
#' @param pii Vector con las probabilidades de inclusión de primer orden para cada elemento de la muestra. Debe tener la misma longitud que \code{y}.
#' @details 
#' La función implementa el estimador de Hajek para la media, estimando el total poblacional multiplicando la media estimada de Hájek por el tamaño de la población \code{N}.
#' @return Una lista con los siguientes componentes:
#' \itemize{
#'    \item \code{Estimación de la media}
#'    \item \code{Estimación del total}
#' }
#' @author Javier Miralles Corbacho
#' @export
#' @examples
#' #Ejemplo aplicado a un muestreo de Bernouilli(N=1000,p=0.005), cuyo vector de muestra es (45, 34, 50, 48, 60, 25)
#' y <- c(45, 34, 50, 48, 60, 25) ; n <- length(y)
#' N <- 1000 ; p <- 0.005
#' pii <- rep(p,n)
#' est_Hajek(y,N,pii)

est_Hajek <- function(y,N,pii)
{
  #validaciones de los parámetros sobre:
  #1.- el vector de muestra 'y': ha de ser vector, la longitud ha de ser menor que el tamaño poblacional N y no debe contener NA valores
  if(is.vector(y)==FALSE) stop("Error, el parámetro 'y' ha de ser un vector.")
  if(length(y)>=N) stop("Error, el tamaño muestral ha de ser inferior al tamaño poblacional 'N'.")
  if(any(is.na(y))==TRUE) stop("Error, el vector de muestra 'y' contiene NA valores.")
  
  #2.- el tamaño poblacional 'N': ha de ser un número natural mayor que 2
  if(N<2) stop("Error, el tamaño poblacional ha de ser un número natural mayor que 2.")

  #3.- el vector de probabilidades de inclusión de primer orden 'pii': ha de ser un vector, estrictamente positivo, debe tener igual longitud que el vector de muestra y deben ser valores inferiores estrictamente que 1
  if(is.vector(pii)==FALSE) stop("Error, el parámetro 'pii' ha de ser un vector.")
  if(any(pii<=0)) stop("Error, las probabilidades de inclusión deben ser estrictamente positivas.")
  if(length(y) != length(pii)) stop("Error, la longitud de la muestra 'y' y del vector de probabilidades 'pii' deben coincidir.")
  if(any(pii>=1)) stop("Error, las probabilidades de inclusión no pueden ser mayores o iguales que 1.")
  
  #calculamos ambos estimadores, previamente calculando los pesos
  w <- 1/pii
  estMedia <- sum(y*w)/sum(w)
  estTotal <- N*estMedia
  
  #devolvemos los valores
  return(list("Estimación de la media"=estMedia,"Estimación del total"=estTotal))
}