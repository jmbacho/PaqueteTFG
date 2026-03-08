#' @title Estimador de Hansen-Hurwitz para la media y el total
#' @description Calcula el estimador de Hansen-Hurwitz del total y de la media poblacional, así como sus varianzas estimadas e intervalos de confianza, generalmente utilizado en diseños con reemplazamiento.
#' @param y Vector con los valores de la muestra observada. No debe contener NAs.
#' @param N Valor que indica el tamaño de la población.
#' @param pi Vector de dimensión n (tamaño muestral) con las probabilidades de selección (o extracción) de cada unidad observada.
#' @param conf_level Nivel de confianza para los intervalos (por defecto 0.95).
#' @details 
#' La función implementa el estimador insesgado de Hansen-Hurwitz para el total y la media..
#' @return Una lista con los siguientes componentes:
#' \itemize{
#'    \item \code{Estimación de la media}
#'    \item \code{Estimación de la varianza media}
#'    \item \code{Intervalo de confianza para la media}
#'    \item \code{Estimación del total}
#'    \item \code{Estimación de la varianza total}
#'    \item \code{Intervalo de confianza para el total}
#' }
#' @author Javier Miralles Corbacho
#' @export
#' @examples
#' #Estimar la producción de leche total y leche media en una aldea con 4 granjas, seleccionando las granjas con probabilidad proporcional al número de vacas; Granja A: 10 vacas, Granja B: 20 vacas, Granja C: 30 vacas, Granja D: 40 vacas.
#' #Tamaño poblacional
#' N <- 4
#' #Probabilidades de selección
#' pi <- c(0.10,0.20,0.30,0.40)
#' #Tomamos una muestra aleatoria simple de tamaño 2, resultando las granjas B y D
#' #Se observan que en las granjas seleccionadas se han recogido 50 y 180 litros de leche respectivamente
#' y <- c(50,180)
#' pi <- pi[c(2,4)]
#' est_HH(y,N,pi)

est_HH <- function(y,N,pi,conf_level=0.95)
{
  #validaciones de los parámetros sobre:
  #1.- el vector de muestra 'y': ha de ser un vector y no puede contener NA valores.
  if(is.vector(y)==FALSE) stop("Error, el parámetro 'y' ha de ser un vector.")
  if(any(is.na(y))==TRUE) stop("Error, el vector de muestra 'y' contiene NA valores.")
  
  #2.- el tamaño poblacional 'N': que sea mayor que el tamaño muestral y que sea un valor mayor o igual que 2
  if((N<1)==TRUE) stop("Error, el tamaño poblacional ha de ser un número positivo mayor que 2.")
  if(length(y)>N) stop("Error, el tamaño muestral ha de ser inferior al tamaño poblacional 'N'.")
  
  #3.- el vector de probabilidades de selección 'pi': que sea un vector, que tenga igual longitud que el vector muestral y que los valores sean mayores estrictos que cero
  if(is.vector(pi)==FALSE) stop("Error, el parámetro 'pi' ha de ser un vector.")
  if(length(y)!=length(pi)) stop("Error, la longitud del vector de probabilidades 'pi' debe coincidir con la del vector muestral 'y'.")
  if(any(pi<=0) || any(pi>1)) stop("Error, las probabilidades de selección deben estar entre 0 y 1.")
  
  #4.- el nivel de confianza 'conf_level': ha de estar entre 0 y 1
  if(conf_level<=0 || conf_level>=1) stop("Error, el nivel de confianza debe ser un valor del intervalo abierto (0,1).")
  
  #determinamos la longitud muestral
  n <- length(y)
  #realizamos los cocientes y_i/p_i (cambio de variable)
  z <- y/pi
  
  #calculamos ambos estimadores
  estTotal <- sum(z)/n
  estMedia <- estTotal/N
  
  #calculamos la estimación de la varianza total
  varEstTotal <- 1/(n*(n-1)) * sum((z - estTotal)^2)
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
