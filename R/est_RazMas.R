#' @title Estimador de razón en Muestreo Aleatorio Simple (MAS)
#' @description Calcula las estimaciones de razón, media y total, así como sus varianzas estimadas utilizando el estimador de razón. Esta función es específica para diseños de Muestreo Aleatorio Simple sin reemplazamiento, utilizando la correlación de la variable de interés con una variable auxiliar.
#' @param y Vector muestral de la variable de interés.
#' @param x Vector muestral de la variable auxiliar. Debe tener la misma longitud que \code{y}.
#' @param tX Valor que indica el total poblacional conocido de la variable auxiliar X.
#' @param N Valor que indica el tamaño total de la población.
#' @details
#' La función implementa el estimador de razón teniendo en cuenta que es sesgado, por lo que también devuelve una estimación del sesgo para evaluar si es o no despreciable.
#' @return Una lista con los siguientes componentes:
#' \itemize{
#'   \item \code{Estimador de razón}
#'   \item \code{Estimación del sesgo del estimador de razón}
#'   \item \code{Estimación de la media}
#'   \item \code{Estimación de la varianza media}
#'   \item \code{Estimación del total}
#'   \item \code{Estimación de la varianza total}
#' }
#' @author Javier Miralles Corbacho
#' @export
#' @examples
#' #datos poblacionales
#' N <- 750   #tamaño poblacional (socios)
#' tX <- 3840 #total poblacional de X (hectáreas)
#' #datos muestrales
#' y <- c(12,14,11,15,16,12,24,15,18,20,8,20,16,14,18,15,18,17,15,22)        #producción
#' x <- c(3.7,4.3,4.1,5,5.5,3.8,8,5.1,5.7,6,3,7,5.4,4.4,5.5,5,5.9,5.6,5,7.2) #superficie
#' est_RazMas(y,x,tX,N)

est_RazMas <- function(y,x,tX,N)
{
  #validaciones de los parámetros sobre:
  #1.- el vector de muestra 'y': ha de ser vector, la longitud ha de ser menor que el tamaño poblacional N y no debe contener NA valores
  if(is.vector(y)==FALSE) stop("Error, el parámetro 'y' ha de ser un vector.")
  if(length(y)>=N) stop("Error, el tamaño muestral ha de ser inferior al tamaño poblacional 'N'.")
  if(any(is.na(y))==TRUE) stop("Error, el vector de muestra 'y' contiene NA valores.")
  
  #2.- el vector de muestra auxiliar 'x': ha de ser vector, no debe contener NAs y debe tener la misma longitud que 'y'
  if(is.vector(x)==FALSE) stop("Error, el parámetro 'x' ha de ser un vector.")
  if(any(is.na(x))==TRUE) stop("Error, el vector de muestra 'x' contiene NA valores.")
  if(length(x)!=length(y)) stop("Error, los vectores muestrales 'y' y 'x' han de tener la misma longitud.")
  
  #3.- el tamaño poblacional 'N': ha de ser un número natural mayor que 2
  if(N<2) stop("Error, el tamaño poblacional ha de ser un número natural mayor que 2.")
  
  #4.- el total poblacional 'tX': ha de ser un número
  if(is.numeric(tX)==FALSE) stop("Error, el parámetro 'tX' ha de ser un valor numérico.")
  
  #determinamos el tamaño muestral
  n <- length(y)
  #determinamos el cociente del tamaño muestral por el tamaño poblacional
  f <- n/N
  
  #calculamos la estimación de la razón R
  estR <- mean(y)/mean(x)
  #calculamos ahora ambas estimaciones
  estTotal <- estR*tX
  estMedia <- estTotal/N
  
  #calculamos las estimaciones de la varianza correspondientes
  varEstMedia <- (1-f)/n * (var(y)+estR^2*var(x)-2*estR*cov(x,y))
  varEstTotal <- varEstMedia*N^2
  
  #sabemos que el estimador de razon es sesgado, por lo que estimemos este sesgo
  estSesgoR <- (1-f)/(n*((tX/N)^2)) * (estR*var(x)-cov(x,y))
  
  #devolvemos los valores
  return(list("Estimador de razón"=estR,
              "Estimación del sesgo del estimador de razón"=estSesgoR,
              "Estimación de la media"=estMedia,
              "Estimación de la varianza media"=varEstMedia,
              "Estimación del total"=estTotal,
              "Estimación de la varianza total"=varEstTotal))
}