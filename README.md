## llrRcpp


### Installation
Install devtools 

```
library(devtools)
devtools::install_github("Chengwei94/llrRcpp")
```

### Details 
Making use of a kd-tree for estimation of local linear estimation. There is an exact kd-tree and an approximate method. Added a metaheuristic approach
to search for the best cv.

### Functions
1. Binning of data 
```
## Binning of 1D/2D data into grids 
x <- runif(n, 0, 1) 
xvar <- rnorm(n, 0, 0.5)
y <- sin(2*pi*x) + xvar
binned <- bin(x, y, bins=400, w)
## output x, y, and weights for the binned data. 
```
2. Local linear regression   
```
## local linear regression for exact
llr_exact <- llr(x, y, x, bandwidth =0.2, weight = w)
## local linear regression for kdtree exact
llr_kdexact <- llr(x, y, x, bandwidth = 0.2, weight = w, kdtree = TRUE)
## local linear regression for kdtree approximation
llr_kdapprox <- llr(x, y, x, bandwidth = 0.2, weight = w, kdtree = TRUE, approx = TRUE)
## local linear regression for data after binning.
llr_bin <- llr(binned, x , bandwidth = 0.2)
```
3. Leave one out Cross validation 
```
## Bandwidth selection of binned data
h_bin <- loocv.llr(binned$x, binned$y, binned$weight)
## Bandwidth selection of exact local linear regression
h_exact <- loocv.llr(x, y, w)
## Bandwidth selection of approx local linear regression with kdtree
h_kdapprox <- loocv.llr(x, y, w, approx = TRUE)
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

5. Autoloocv selection 

```
## Bandwidth selection using auto exact
h_auto <- autoloocv.llr(x, y, w)

## Bandwdith 
h_auto <- autoloocv.llr(x, y, w, approx = TRUE) 
```
### Visual Results (Red stand for bandwidth chosen with cv (approx), Yellow stand for bandwidth chosen gcv(exact), Blue stand for bandwidth chosen with gcv(approx)).
x1 ~ N(0,0.3)  
x ~ U(0,10)  
y = sin(x) + x1  
N = 100  
![1](https://user-images.githubusercontent.com/61018420/103854614-df5eca00-50eb-11eb-9975-713c1f8b6345.jpg)  
N = 500  
![2](https://user-images.githubusercontent.com/61018420/103856245-4c279380-50ef-11eb-8818-43dd68ee6397.jpg)  
N = 1000  
![3](https://user-images.githubusercontent.com/61018420/103856247-4cc02a00-50ef-11eb-873b-c4a9d20f9f58.jpg)  
N = 5000  
![4](https://user-images.githubusercontent.com/61018420/103856248-4cc02a00-50ef-11eb-84a9-e701b283c7ac.jpg)  

### Accuracy Results
1D comparison  
![5](https://github.com/Chengwei94/llrRcpp/blob/main/inst/1d.png)

3D comparison  
![6](https://github.com/Chengwei94/llrRcpp/blob/main/inst/3d.png)

4D comparison  
![7](https://github.com/Chengwei94/llrRcpp/blob/main/inst/4d.png)

