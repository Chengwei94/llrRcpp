# Locallinearregression


# Installation
devtools::install_github("Chengwei94/llrRcpp")

# Functions
1. Binning of data
```
## Binning of 1D/2D data into grids 
binned <- bin(x,y,bins=400, w)
## output x, y, and weights for the binned data. 
```
2. Local linear regression 
```
## local linear regression for exact
llr_exact <- llr(x, y, x, bw =0.2, weight = w)
## local linear regression for kdtree exact
llr_kdexact <- llr(x, y, x, bw = 0.2, weight = w, kdtree = TRUE)
## local linear regression for kdtree approximation
llr_kdapprox <- llr(x, y, x, bw = 0.2, weight = w, kdtree = TRUE, approx = TRUE)
## local linear regression for data after binning.
llr_bin <- llr(binned, x , bw = 0.2)
```
3. Generalized Cross validation 
```
## Bandwidth selection of binned data
h_bin <- gcv.llr(binned$x, binned$y, binned$weight)
## Bandwidth selection of exact local linear regression
h_exact <- gcv.llr(x, y, w)
## Bandwidth selection of approx local linear regression with kdtree
h_kdapprox <- gcv.llr(x, y, w, approx = TRUE)
```

4. K-fold Cross validation 
```
binned <- bin(x,y,bins=400, w)
## Bandwidth selection of binned data
h_bin <- cv.llr(binned$x, binned$y, binned$weight)
## Bandwidth selection of exact local linear regression
h_exact <- cv.llr(x, y, w)
## Bandwidth selection of exact local linear regression with kdtree
h_kdexact <- cv.llr(x, y, w, kdtree = TRUE, approx = FALSE)
## Bandwidth selection of approx local linear regression with kdtree
h_kdapprox <- cv.llr(x, y, w , kdtree = TRUE , approx = TRUE)
```


