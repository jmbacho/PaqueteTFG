#' @title Estimador de Horvitz-Thompson para la media y el total
#' @description Calcula el estimador de Horvitz-Thompson del total y de la media poblacional, así como sus varianzas estimadas e intervalos de confianza, generalmente utilizado en diseños sin reemplazamiento.
#' @param y Vector con los valores de la muestra observada. No debe contener NAs.
#' @param N Valor que indica el tamaño de la población.
#' @param P Matriz cuadrada de dimensión n (tamaño muestral) de probabilidades de inclusión. La diagonal representa las probabilidades de primer orden el resto las de segundo orden.
#' @param conf_level Nivel de confianza para los intervalos (por defecto 0.95).
#' @details 
#' La función implementa el estimador insesgado de Horvitz-Thompson para el total y la media.
#' Para la estimación de la varianza, se utiliza la fórmula clásica de Horvitz-Thompson, por lo que requiere que las probabilidades de inclusión sean estrictamente positivas.
#' @return Una lista con los siguientes componentes:
#' \itemize{
#'   \item \code{Estimación de la media}
#'   \item \code{Estimación de la varianza media}
#'   \item \code{Intervalo de confianza para la media}
#'   \item \code{Estimación del total}
#'   \item \code{Estimación de la varianza total}
#'   \item \code{Intervalo de confianza para el total}
#' }
#' @author Javier Miralles Corbacho
#' @export
#' @examples
#' #Ejemplo aplicado a un muestreo de Bernouilli(N=1000,p=0.005), cuyo vector de muestra es (45, 34, 50, 48, 60, 25)
#' y <- c(45, 34, 50, 48, 60, 25)
#' n <- length(y)
#' N <- 1000
#' p <- 0.005
#' P <- matrix(data=p^2, ncol=n, nrow=n) ; diag(P) <- p
#' est_HT(y,N,P)

est_HT <- function(y,N,P,conf_level=0.95)
{
  #validaciones de los parámetros sobre:
  #1.- el vector de muestra 'y': ha de ser vector, la longitud ha de ser menor que el tamaño poblacional N y no debe contener NA valores
  if(is.vector(y)==FALSE) stop("Error, el parámetro 'y' ha de ser un vector.")
  if(length(y)>=N) stop("Error, el tamaño muestral ha de ser inferior al tamaño poblacional 'N'.")
  if(any(is.na(y))==TRUE) stop("Error, el vector de muestra 'y' contiene NA valores.")
  
  #2.- el tamaño poblacional 'N': ha de ser un número natural mayor que 2
  if(N<2) stop("Error, el tamaño poblacional ha de ser un número natural mayor que 2.")
  
  #3.- la matriz de probabilidades de inclusión 'P': que sea una matriz, que todas los valores sean mayores estrictos que cero y que tenga las dimensiones adecuadas
  if(is.matrix(P)==FALSE) stop("Error, el parámetro P ha de ser una matriz.")
  if(any(P<=0)) stop("Error, el diseño no es cuantificable.")
  if(ncol(P)!=nrow(P) ||  nrow(P)!=length(y)) stop("Error, la matriz P debe ser cuadrada y tener dimensión igual al tamaño muestral n.")
  
  #4.- el nivel de confianza 'conf_level': ha de estar entre 0 y 1
  if(conf_level<=0 || conf_level>=1) stop("Error, el nivel de confianza debe ser un valor del intervalo abierto (0,1).")
  
  #tamaño muestral
  n <- length(y)
  
  #extraemos de la matriz de probabilidades P las probabilidades de inclusión de primer orden
  pii <- diag(P)
  
  #calculamos las estimaciones de Horvitz-Thompson del total y de la media
  estTotal <- sum(y/pii)
  estMedia <- estTotal/N
  
  #calculamos la estimación de la varianza total (insesgada), para ello primero calculamos los delta_ij
  Dij <- P - outer(pii,pii)
  #ahora calculamos los cocientes y_i/pi_i
  yp <- y/pii
  #calculamos los cocientes delta_ij/pi_ij
  dp <- Dij/P
  #realizamos el producto matricial utilizando %*%
  varEstTotal <- as.numeric(t(yp) %*% dp %*% yp)
  if(varEstTotal<0)
  {
    warning("El cálculo de la varianza estimada para el total es negativa, por lo que se devuelve el valor cero.")
    varEstTotal=0
  }
  
  #calculamos ahora la estimación de la varianza media (insesgada)
  varEstMedia <- varEstTotal/N^2
  
  #calculamos intervalos de confianza para la estimación del total y de la media
  zt <- qt(1-(1-conf_level)/2, n-1)*sqrt(varEstTotal)
  it <- c(inf=estTotal-zt, sup=estTotal+zt)
  zm <- qt(1-(1-conf_level)/2, n-1)*sqrt(varEstMedia)
  im <- c(inf=estMedia-zm, sup=estMedia+zm)
  
  #creamos variables dinámicas para los intervalos de confianza
  ci_media <- paste0("Intervalo de confianza al nivel ", conf_level, " para la media")
  ci_total <- paste0("Intervalo de confianza al nivel ", conf_level, " para el total")
  
  #devolvemos los valores
  return(
    setNames(
      list(estMedia, varEstMedia, im, estTotal, varEstTotal, it),
      c("Estimación de la media",
        "Estimación de la varianza media",
        ci_media,
        "Estimación del total",
        "Estimación de la varianza total",
        ci_total)
    )
  )
}
