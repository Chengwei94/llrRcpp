#' Plot class 'llr' object
#' @rdname plot
#' @method plot llr 
#' @param x a object of class 'llr' 
#' @param xorig initial numeric vector x 
#' @param yorig initial numeric vector y
#' @param ... for consistency
#' @export
plot.llr <- function(x, xorig, yorig,  ...) {
  plot(xorig, yorig) 
  lines(x$x, x$fitted, col = 'red')
}

#' Plot 'bin' object
#' @rdname plot 
#' @method plot bin 
#' @param x a object of class 'bin'
#' @param ... for consistency 
#' @export
plot.bin <- function(x, ...) {
  plot(x$x, x$y)
}