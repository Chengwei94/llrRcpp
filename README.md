## llrRcpp


## Installation
Install devtools 

```
library(devtools)
devtools::install_github("Chengwei94/llrRcpp")
```

# Details 
Making use of a kd-tree for estimation of local linear estimation. For both approx and exact method. 

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
# Visual Results (Red stand for bandwidth chosen with cv (approx), Yellow stand for bandwidth chosen gcv(exact), Blue stand for bandwidth chosen with gcv(approx)).
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

# Microbenchmark Speed Results (With similar bandwidth)  
1. LLR speed (compared to Kernsmooth Locpoly(only for 1D) and loess)   
A. 1D  
y = sin(x) + x1  
x1 ~ N(0, 0.3)  
x ~ U(0, 10)  
N = 100  
![mb-0](https://user-images.githubusercontent.com/61018420/103856255-4df15700-50ef-11eb-8374-541801c7bd6b.JPG)  
N = 1000  
![mb-1](https://user-images.githubusercontent.com/61018420/103856233-492ca300-50ef-11eb-8899-5763b763ebfa.JPG)  
N = 10000  
![mb-3](https://user-images.githubusercontent.com/61018420/103856237-4a5dd000-50ef-11eb-8edb-9986f0763628.JPG)  
  

# B. 2D  
N = 1000  
![mb-4](https://user-images.githubusercontent.com/61018420/103856239-4af66680-50ef-11eb-8035-0206a9f1b9eb.JPG)   
N = 10000  
![mb-5](https://user-images.githubusercontent.com/61018420/103856240-4b8efd00-50ef-11eb-8213-1bda98f8b886.JPG)   


# C. 3D   
N = 1000   
![bm-6](https://user-images.githubusercontent.com/61018420/103856250-4d58c080-50ef-11eb-85ea-2546f1997c39.JPG)  
N = 10000  
![bm-7](https://user-images.githubusercontent.com/61018420/103856252-4df15700-50ef-11eb-9dd7-07630b62c623.JPG)  

# 2. Cross-validation time  
N = 1000    
![cv-1](https://user-images.githubusercontent.com/61018420/103857111-d15f7800-50f0-11eb-96db-6551d19c65ac.JPG)  
N = 5000  
![cv-2](https://user-images.githubusercontent.com/61018420/103857164-e89e6580-50f0-11eb-8bd9-d963bfc5b49a.JPG)  

# Accuracy Results
To be moved over 

