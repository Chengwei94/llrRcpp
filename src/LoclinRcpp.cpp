#include "RcppEigen.h"
#include "LoclinRcpp_types.h"
#include "LoclinRcpp.h"
#include <utility>
#include <algorithm>
#include <iostream>   
#include <random>
#include <stack>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
  
// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

//Constructor for kdtree
kdtree::kdtree(const Eigen::MatrixXd& X, const Eigen::VectorXd& Y, const Eigen::VectorXd& wt,  int N_min) {   // command to create kd tree 
  
  Eigen::MatrixXd XY_mat(X.rows(), X.cols() + 2);
  XY_mat << X, Y, wt; 
  all_point_t XY_arr = XYmat_to_XYarr(XY_mat) ;  
  
  size_t len = XY_arr.size(); 
  
  std::vector<double> dim_max_(X.cols());        
  std::vector<double> dim_min_(X.cols()); 
  
  for(size_t i=0; i < X.cols(); i++) { 
    dim_max_.at(i) = X.col(i).maxCoeff(); 
    dim_min_.at(i) = X.col(i).minCoeff(); 
  } 
  
  root = build_tree(XY_arr.begin(), XY_arr.end(), 1, 1, N_min, len, dim_max_, dim_min_);
}

// Function for building tree, storing max_dim , min_dim, split dimension, split point, length at each node.
std::unique_ptr<kdtree::kdnode> kdtree::build_tree(all_point_t::iterator start, all_point_t::iterator end, 
                                           int split_d, double split_v, int N_min, size_t len,
                                           std::vector<double> dim_max_, std::vector<double> dim_min_){
  
  std::unique_ptr<kdtree::kdnode> newnode = std::make_unique<kdtree::kdnode>();
  
  newnode->n_below = len;    
  
  if (end-start <= N_min) {            
    Eigen::MatrixXd XtX_(dim_max_.size() + 1 , dim_max_.size() + 1);
    XtX_.setZero(); 
    Eigen::VectorXd XtY_(dim_max_.size() + 1);
    XtY_.setZero();
    for (size_t k = 0; k <dim_max_.size(); k++ ){
      dim_max_[k] = (*start)(k+1);
      dim_min_[k] = (*start)(k+1);
    }                  
    for (auto i = start; i != end; i ++){
      Eigen::VectorXd XYw = *i;   
      double wt = XYw(XYw.size()-1);
      double Y = XYw(XYw.size()-2);
      Eigen::VectorXd X = XYw.head(XYw.size()-2); 
      XtY_ += X * Y * wt; 
      XtX_ += wt * X * X.transpose();  
      for (size_t j = 0; j < dim_max_.size(); j++) {
        if(X(j+1) > dim_max_[j]){
          dim_max_[j] = X(j+1); 
        }
        if(X(j+1) < dim_min_[j]){
          dim_min_[j] = X(j+1); 
        } 
      }
    } 
    newnode->dim_max = dim_max_;
    newnode->dim_min = dim_min_;   
    newnode->XtX = XtX_; 
    newnode->XtY = XtY_;  
    return newnode; 
  }
  
  else {
    newnode->dim_max = dim_max_;
    newnode->dim_min = dim_min_;  
    
    size_t l_len = len/2  ;          
    size_t r_len = len - l_len;    
    auto middle = start + len/2;   
    int max = 0; 
    int dim = 0;
    for(size_t i = 0; i < newnode->dim_max.size(); i++){   
      double var = newnode->dim_max[i] - newnode->dim_min[i]; 
      if(var > max){
        max = var; 
        dim = i; 
      }
    }
    
    newnode->split_d = dim; 
    int vector_dim = dim + 1;  
    
    std::nth_element(start, middle, end, [vector_dim](const Eigen::VectorXd& a, const Eigen::VectorXd & b) {
      return a(vector_dim) < b(vector_dim);    
    });           
    
    newnode->split_v = (*middle)[vector_dim];   
    dim_max_[dim] = newnode->split_v; 
    dim_min_[dim] = newnode->split_v;
    
    newnode-> left_child = build_tree(start, middle, newnode->split_d, newnode->split_v, N_min, l_len, dim_max_, newnode->dim_min);
    newnode-> right_child = build_tree(middle, end, newnode->split_d, newnode->split_v, N_min, r_len, newnode->dim_max, dim_min_);
    
    if ((newnode->left_child) && (newnode->right_child)){ 
      newnode->XtX = newnode->left_child->XtX + newnode ->right_child->XtX;  // sumY = the sum of the bottom 2 nodes
      newnode->XtY = newnode->left_child->XtY + newnode ->right_child->XtY;
    } 
    else if (newnode->left_child) {
      newnode->XtY = newnode->left_child->XtY;
      newnode->XtX = newnode->left_child->XtX;
    }
    else if (newnode->right_child) {
      newnode->XtX = newnode->right_child->XtX; 
      newnode->XtY = newnode->right_child->XtY; 
    }   
  }
  return newnode;    
}

// Function to change from matrix form to std::vector form since that row-wise sorting is not 
// supported in Eigen
all_point_t XYmat_to_XYarr(const Eigen::MatrixXd& XY_mat){    
  std::vector<Eigen::VectorXd> XY_arr; 
  Eigen::MatrixXd XY_temp(XY_mat.rows(), XY_mat.cols()+1); 
  XY_temp.rightCols(XY_mat.cols()) = XY_mat; 
  XY_temp.col(0) = Eigen::VectorXd::Ones(XY_mat.rows()); 
  
  for (int i = 0; i <XY_temp.rows(); i++){
    XY_arr.push_back(XY_temp.row(i));
  }
  return XY_arr; 
}

// evaluation of kernel
inline double eval_kernel(int kcode, const double& z, double w)
{
  double tmp;
  if(abs(z) > 1) 
    return 0;
  else { 
    switch(kcode)
    {
    case 1: return (3*(1-z*z)/4)*w; // Epanechnikov
    case 2: return 0.5*w; // rectangular
    case 3: return (1-abs(z)) * w; // triangular
    case 4: return (15*(1-z*z)*(1-z*z)/16) * w; // quartic
    case 5: 
      tmp = 1-z*z;
      return (35*tmp*tmp*tmp/32)*w; // triweight
    case 6: 
      tmp = 1- z*z*z;
      return (70 * tmp * tmp * tmp / 81)*w; // tricube
    case 7:
      return (M_PI * cos(M_PI*z/2) / 4)*w; // cosine
    case 21:
      return (exp(-z*z/2) / sqrt(2*M_PI))*w; // gauss
    case 22:
      return (1/(exp(z)+2+exp(-z)))*w; // logistic
    case 23:
      return (2/(M_PI*(exp(z)+exp(-z))))*w; // sigmoid
    case 24:
      return (exp(-abs(z)/sqrt(2)) * sin(abs(z)/sqrt(2)+M_PI/4))*w; // silverman
      //   default: Rcpp::stop("Unsupported kernel"); 
    }
  }
  return 0;
}

inline double eval_kernel(int kcode, const double& z)
{
  double tmp;
  if(abs(z) > 1) 
    return 0;
  else { 
    switch(kcode)
    {
    case 1: return 3*(1-z*z)/4; // Epanechnikov
    case 2: return 0.5; // rectangular
    case 3: return 1-abs(z); // triangular
    case 4: return 15*(1-z*z)*(1-z*z)/16; // quartic
    case 5: 
      tmp = 1-z*z;
      return 35*tmp*tmp*tmp/32; // triweight
    case 6: 
      tmp = 1- z*z*z;
      return 70 * tmp * tmp * tmp / 81; // tricube
    case 7:
      return M_PI * cos(M_PI*z/2) / 4; // cosine
    case 21:
      return exp(-z*z/2) / sqrt(2*M_PI); // gauss
    case 22:
      return 1/(exp(z)+2+exp(-z)); // logistic
    case 23:
      return 2/(M_PI*(exp(z)+exp(-z))); // sigmoid
    case 24:
      return exp(-abs(z)/sqrt(2)) * sin(abs(z)/sqrt(2)+M_PI/4); // silverman
      //   default: Rcpp::stop("Unsupported kernel"); 
    }
  }
  return 0;
}

// calculate max weight and min weight between the query point and its bounding box
std::pair<double,double> calculate_weight(int kcode, const Eigen::VectorXd& X_query, const std::vector<double>& dim_max, 
                                          const std::vector<double>& dim_min, const Eigen::VectorXd& h) {  // weight calculation of point
  if(dim_max.size()!= X_query.size()){
    //  std::cout << "error in dimensions"; 
    //    throw std::exception(); 
  }
  double w_max = 1;
  double w_min = 1;
  
  for(int i=0; i < dim_max.size(); i++) { 
    if (X_query(i) <= dim_max.at(i) && X_query(i) >= dim_min.at(i)) {
      double dis_max = std::max(abs(dim_max.at(i) - X_query(i)), abs(dim_min.at(i) - X_query(i)));    // dis_max = which ever dim is further
      double dis_min = 0;                                                                             
      w_min *= eval_kernel(kcode, dis_max/h(i)) / h(i);    
      w_max *= eval_kernel(kcode, dis_min/h(i)) / h(i);
    }
    else{
      double dist1 = abs(X_query(i) - dim_min.at(i));
      double dist2 = abs(X_query(i) - dim_max.at(i));
      double dis_min = std::min(dist1, dist2);
      double dis_max = std::max(dist1, dist2);
      w_min *= eval_kernel(kcode, dis_max/h(i)) / h(i); 
      w_max *= eval_kernel(kcode, dis_min/h(i)) / h(i); 
    }
  }
  return std::make_pair(w_max, w_min); 
}

// Function to get exact XtX XtY based on condition on w_max and w_min is equal
std::pair<Eigen::MatrixXd, Eigen::VectorXd> kdtree::get_XtXXtY(const Eigen::VectorXd& X_query,  // to obtain the exact XtXXtY from tree 
                                                               const std::vector<double>& dim_max, 
                                                               const std::vector<double>& dim_min, 
                                                               std::unique_ptr<kdtree::kdnode>& root, 
                                                               const Eigen::VectorXd& h, 
                                                               int kcode){
  
  std::pair<double,double> w_maxmin;
  std::stack<kdtree::kdnode*> storage; 
  auto curr = root.get(); 
  Eigen::MatrixXd XtX = Eigen::MatrixXd::Zero(curr->XtX.rows() , curr->XtX.cols()); 
  Eigen::VectorXd XtY = Eigen::MatrixXd::Zero(curr->XtY.rows() , curr->XtY.cols());
  w_maxmin = calculate_weight(kcode, X_query, dim_max, dim_min, h);
  double w_max = w_maxmin.first; 
  double w_min = w_maxmin.second;
  double count; 
  
  while (w_max != w_min || storage.empty() == false){
    while (w_max != w_min ){   
      count += 1; 
      storage.push(curr);
      curr = curr->left_child.get();  
      w_maxmin = calculate_weight(kcode, X_query, curr->dim_max, curr->dim_min, h);   
      w_max = w_maxmin.first;                      
      w_min = w_maxmin.second; 
      
      if(w_max == w_min ){ 
        XtX += w_max*curr->XtX;
        XtY += w_max*curr->XtY; 
      }           
    }
    count += 1; 
    curr = storage.top();
    storage.pop(); 
    curr = curr->right_child.get(); 
    w_maxmin = calculate_weight(kcode, X_query, curr->dim_max, curr->dim_min, h);  // calculate max and min weight 
    w_max = w_maxmin.first;   
    w_min = w_maxmin.second; 
    if(w_max == w_min){ 
      XtX += w_max*curr->XtX;
      XtY += w_max*curr->XtY;
    }
  }   
  // Rcpp::Rcout << count <<'\n';
  return std::make_pair(XtX,XtY); 
}

//Function to get approx XtX XtY from the tree choosing based on the condition w_max - w_min <= 2 * current_weight + w_min * nodes below
std::pair<Eigen::MatrixXd, Eigen::VectorXd> kdtree::getapprox_XtXXtY(const Eigen::VectorXd& X_query,    
                                                                     const std::vector<double>& dim_max,
                                                                     const std::vector<double>& dim_min, 
                                                                     std::unique_ptr<kdtree::kdnode>& root,
                                                                     double epsilon, 
                                                                     const Eigen::VectorXd& h,
                                                                     int kcode){ 
  
  std::pair<double,double> w_maxmin;
  std::stack<kdtree::kdnode*> storage; 
  auto curr = root.get(); 
  Eigen::MatrixXd XtX = Eigen::MatrixXd::Zero(curr->XtX.rows() , curr->XtX.cols()); 
  Eigen::VectorXd XtY = Eigen::MatrixXd::Zero(curr->XtY.rows() , curr->XtY.cols());
  
  w_maxmin = calculate_weight(kcode, X_query, dim_max, dim_min, h);
  double w_max = w_maxmin.first; 
  double w_min = w_maxmin.second;
  
  double weight_sf = 0; 
  double count = 0;
  
  while ((w_max-w_min > 2*epsilon*(weight_sf + curr->n_below*w_min)) || storage.empty() == false  ){
    while (w_max - w_min > 2*epsilon*(weight_sf + curr->n_below*w_min)){   // if condition fufilled        
      count += 1;
      storage.push(curr);
      curr = curr->left_child.get();  
      w_maxmin = calculate_weight(kcode, X_query, curr->dim_max, curr->dim_min, h);  // calculate max and min weight 
      w_max = w_maxmin.first;                      
      w_min = w_maxmin.second;
      
      if(w_max - w_min <= 2 * epsilon * (weight_sf + curr->n_below*w_min)){ 
        double weight = 0.5*(w_max + w_min);
        //weight_sf += weight;  
        XtX += weight * curr->XtX;
        XtY += weight * curr->XtY; 
      }           
    }
    
    count += 1;
    curr = storage.top(); 
    storage.pop(); 
    curr = curr->right_child.get(); 
    w_maxmin = calculate_weight(kcode, X_query, curr->dim_max, curr->dim_min, h);  // calculate max and min weight 
    w_max = w_maxmin.first;   
    w_min = w_maxmin.second;
    
    if(w_max - w_min <= 2 * epsilon * (weight_sf + curr->n_below*w_min)){ 
      double weight = 0.5*(w_max + w_min);
      //weight_sf += weight;  
      XtX += weight * curr->XtX;
      XtY += weight * curr->XtY; 
    }          
  }
  // Rcpp::Rcout << count <<'\n';
  return std::make_pair(XtX,XtY); 
}

// Function to find the XtXXtY from the kdtree
std::pair<Eigen::MatrixXd, Eigen::VectorXd> kdtree::find_XtXXtY(const Eigen::VectorXd& X_query, 
                                                                 int method, 
                                                                double epsilon, 
                                                                const Eigen::VectorXd& h, 
                                                                int kcode ) { 
  if (method == 1) {
    std::pair<Eigen::MatrixXd, Eigen::VectorXd> XtXXtY = get_XtXXtY(X_query, root->dim_max, root->dim_min, root, h, kcode);
    Eigen::MatrixXd ll_XtX = form_ll_XtX(XtXXtY.first, X_query);
    Eigen::VectorXd ll_XtY = form_ll_XtY(XtXXtY.second, X_query);
    return std::make_pair(ll_XtX , ll_XtY); 
  }
  else {
    std::pair<Eigen::MatrixXd, Eigen::VectorXd> XtXXtY = getapprox_XtXXtY(X_query, root->dim_max, root->dim_min, root, epsilon, h, kcode); 
    Eigen::MatrixXd ll_XtX = form_ll_XtX(XtXXtY.first, X_query); 
    Eigen::VectorXd ll_XtY = form_ll_XtY(XtXXtY.second, X_query);
    return std::make_pair(ll_XtX , ll_XtY); 
  }
}

// Function to form XtX into the one which we can use for local linear regression
Eigen::MatrixXd form_ll_XtX(const Eigen::MatrixXd& XtX, const Eigen::VectorXd& X_query){  //steps to form the local linear XtX(with query pts)
  Eigen::MatrixXd extra_XtX(XtX.rows(), XtX.cols()); 
  Eigen::MatrixXd ll_XtX(XtX.rows(),XtX.cols()); 
  extra_XtX.topLeftCorner(1,1) = Eigen::MatrixXd::Zero(1,1); 
  extra_XtX.topRightCorner(1,XtX.cols()-1) = XtX.topLeftCorner(1,1) * X_query.transpose(); 
  extra_XtX.bottomLeftCorner(XtX.rows()-1, 1) =  X_query * XtX.topLeftCorner(1,1);
  extra_XtX.bottomRightCorner (XtX.rows()-1, XtX.cols()-1) = XtX.bottomLeftCorner(XtX.rows()-1,1)*X_query.transpose() 
    + X_query * XtX.topRightCorner(1,XtX.cols()-1) 
    - X_query * XtX.topLeftCorner(1,1) * X_query.transpose();                                             
    ll_XtX = XtX - extra_XtX; 
    return ll_XtX; 
}

// Function to form XtY into the one which we can use for local linear regression
Eigen::VectorXd form_ll_XtY(const Eigen::VectorXd& XtY, const Eigen::VectorXd& X_query){ //steps to form the local linear XtY(with query pts)
  Eigen::VectorXd extra_XtY = Eigen::VectorXd::Zero((XtY.size()));
  Eigen::VectorXd ll_XtY(XtY.size());
  extra_XtY.tail(XtY.size()-1) = X_query * XtY.head(1); 
  ll_XtY = XtY - extra_XtY; 
  return ll_XtY;
}

double calculate_beta(const Eigen::MatrixXd &XtX, const Eigen::VectorXd &XtY) {
  Eigen::LLT<Eigen::MatrixXd> LLT_XtX(XtX);
  if (LLT_XtX.info() == Eigen::NumericalIssue){ 
     return R_PosInf;
   }
  else{ 
//  Eigen::VectorXd beta = XtX.ldlt().solve(XtY);
    Eigen::VectorXd beta = LLT_XtX.solve(XtY); // returns b = (XTX)^-1(XTY)
    return (beta(0));
   }
}
// Function to find the max weight possible for a given node
double max_weight(int kcode, const Eigen::VectorXd &h, double wt){ //return max_weight possible of a point 
  double w_max = 1; 
  for(size_t i =0; i < h.size(); i++){
    w_max *= eval_kernel(kcode, 0)/h(i);
  }
  return w_max * wt;
}

// [[Rcpp::export]]
Eigen::VectorXd llr1d_cpp(const Eigen::VectorXd& X, const Eigen::VectorXd& Y, const Eigen::VectorXd& X_pred, int kcode, double h, Eigen::VectorXd wt){
  Eigen::VectorXd beta_0(X_pred.size());
  int lb = 0;
  int ub = 0; 
  for(int i= 0; i < X_pred.size(); i++) { 
    while (lb < X.size() &&  X(lb) < (X_pred(i) - h)){
      lb++; 
    } 
    while (ub < X.size() &&  X(ub) < (X_pred(i) + h)){ 
      ub++;
    }
    int m = ub - lb;
    
    double s0 = 0 ;
    double s1 = 0 ;
    double s2 = 0 ;
    double t0 = 0 ;
    double t1 = 0 ;
    
    for (int j = 0; j < m; j++){
       double dif = X(j+lb) - X_pred(i);
       double w = eval_kernel(kcode, dif/h)/h * wt(j+lb);
       s0 += w;
       s1 += w * dif;
       s2 += w * pow(dif,2);
       t0 += w * Y(j+lb);
       t1 += w * dif * Y(j+lb);
    }
     double beta = (s2*t0 - s1*t1)/(s0*s2 - pow(s1,2)); 
     if (isnan(beta)){
      beta = R_PosInf; 
     } 
     beta_0(i) = beta;
  }
  return beta_0; 
}

// [[Rcpp::export]]
Eigen::VectorXd llr2d_cpp (const Eigen::MatrixXd &X, const Eigen::VectorXd &Y, const Eigen::MatrixXd& X_pred, int kcode, Eigen::VectorXd h, Eigen::VectorXd wt){
  Eigen::VectorXd beta_0(X_pred.rows());
  int lb = 0;
  int ub = 0; 
  for(int i= 0; i < X_pred.rows(); i++) { 
    while (lb < X.rows() &&  X(lb,0) < (X_pred(i,0) - h(0))){
      lb++; 
    } 
    while (ub < X.rows() &&  X(ub,0) < (X_pred(i,0) + h(0))){ 
      ub++;
    }
    
    Eigen::MatrixXd s(3,3); 
    s.setZero();
    Eigen::VectorXd t(3); 
    t.setZero();  
  
    for (int j = lb; j < ub; j++) { 
      if (abs(X(j,1) - X_pred(i,1)) <= h(1)){
        double x1 = X(j,0) - X_pred(i,0);
        double x2 = X(j,1) - X_pred(i,1);
        double w = wt(j) * eval_kernel(kcode,x1/h(0))/h(0) * eval_kernel(kcode, x2/h(1))/h(1);
        s(0,0) += w; 
        s(0,1) += w * x1; 
        s(1,0) += w * x1;
        s(0,2) += w * x2;
        s(2,0) += w * x2; 
        s(1,1) += w * pow(x1, 2); 
        s(1,2) += w * x1 * x2; 
        s(2,1) += w * x1 * x2 ;
        s(2,2) += w * pow(x2, 2); 
        t(0) += w * Y(j);
        t(1) += w * x1 * Y(j); 
        t(2) += w * x2 * Y(j);
      }
    }
    beta_0(i) = calculate_beta(s, t);
  }
  return beta_0; 
}

// Use to check the kdexact tree
// [[Rcpp::export]]
Eigen::VectorXd llr_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& Y, const Eigen::MatrixXd& X_pred,
                        int kcode, Eigen::VectorXd h, const Eigen::VectorXd &wt){
  Eigen::VectorXd beta_0(X_pred.rows());
  for (int i=0; i < X_pred.rows(); i++){
    Eigen::MatrixXd s(X.cols()+1, X.cols()+1);
    s.setZero();
    Eigen::VectorXd t(X.cols()+1);
    t.setZero(); 
    for (int j=0; j < X.rows(); j++){
      Eigen::VectorXd Xdif = X.row(j) - X_pred.row(i); 
      double w = wt(j); 
      for (int k = 0; k < X.cols(); k++) {
        w *= eval_kernel(kcode, Xdif(k)/h(k))/h(k); 
      }
      s(0,0) += w; 
      t(0) += w * Y(j);
      for (int m = 1; m < X.cols()+1; m++) {
        s(0,m) += w * Xdif(m-1); 
        t(m) +=  w * Xdif(m-1) * Y(j); 
        for (int n = 1; n <= m; n++) { 
        s(n,m) += w * Xdif(m-1) * Xdif(n-1); 
        }
      }
    } 
    s.triangularView<Eigen::Lower>() = s.transpose();
    beta_0(i) = calculate_beta(s, t);
  }
  return beta_0; 
}

// [[Rcpp::export]]
Rcpp::List bin1d_cpp(const Eigen::VectorXd& X, const Eigen::VectorXd&Y, int bins, const Eigen::VectorXd& wt){
  
  double grid_start = X.minCoeff(); 
  double grid_end = X.maxCoeff(); 
  double binwidth = (grid_end - grid_start) / bins;
  
  //if (X.size() < bins || h < binwidth){ 
  //  Rcpp::Rcout << "bandwidth too small or datapoints > bins, using non-bin local linear instead."; 
  //   return predict1dd(X, Y, X_pred, kcode, h);  
  //}
  
  // Linear binning;
  Eigen::VectorXd c_l = Eigen::VectorXd::Zero(bins+1);  
  Eigen::VectorXd d_l = Eigen::VectorXd::Zero(bins+1); 
  Eigen::VectorXd xgrid = Eigen::VectorXd::Zero(bins+1); 
  for (int i = 0; i < bins + 1; i ++){ 
    xgrid(i) = grid_start + i*binwidth; 
  }
  // add weights; 
  for (int i=0; i < X.size() ; i++) {
    double val = (X(i) - grid_start)/binwidth;
    int val_int = floor(val); 
    double val_rem = val - val_int; 
    c_l(val_int) += (1 - val_rem) * wt(i); 
    d_l(val_int) += (1 - val_rem) * Y(i) * wt(i);
    
    if (val_int < bins) {
      c_l(val_int + 1) += val_rem * wt(i); 
      d_l(val_int + 1) += val_rem * Y(i) * wt(i); 
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("weight") = c_l, 
                            Rcpp::Named("y") = d_l.array()/c_l.array(),
                            Rcpp::Named("x") = xgrid);
}

// [[Rcpp::export]]
Rcpp::List bin2d_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& Y, int bins, const Eigen::VectorXd& wt){ // bins at each dir
  
  double grid_start1 = X.col(0).minCoeff(); 
  double grid_end1 = X.col(0).maxCoeff();
  double grid_start2 = X.col(1).minCoeff(); 
  double grid_end2 = X.col(1).maxCoeff(); 
  double binwidth1 = (grid_end1 - grid_start1) / bins;
  double binwidth2 = (grid_end2 - grid_start2) / bins;
  double binarea = binwidth1 * binwidth2; 
  
  //if (X.size() < bins || h < binwidth){ 
  //  Rcpp::Rcout << "bandwidth too small or datapoints > bins, using non-bin local linear instead."; 
  //   return predict1dd(X, Y, X_pred, kcode, h);  
  //}
  
  // Linear binning;
  double ngrid = (bins+1)*(bins+1);
  Eigen::VectorXd c_l = Eigen::VectorXd::Zero(ngrid);  
  Eigen::VectorXd d_l = Eigen::VectorXd::Zero(ngrid); 
  Eigen::MatrixXd xgrid = Eigen::MatrixXd::Zero(ngrid,2); 
  for (int i = 0; i < bins+1; i++){
    for(int j = 0; j < bins+1; j++){
      xgrid(i*(bins+1) + j, 0) = grid_start1 + i*binwidth1; 
      xgrid(i*(bins+1) + j, 1) = grid_start2 + j*binwidth2; 
    }
  }
  // add weights; 
  for (int i=0; i < X.rows(); i++) {
    double val1 = (X(i,0) - grid_start1)/binwidth1;
    double val2 = (X(i,1) - grid_start2)/binwidth2; 
  
    int val_int1 = floor(val1);
    int val_int2 = floor(val2);
    
    double val_rem1 = val1 - val_int1;
    double val_rem2 = val2 - val_int2; 
    
    c_l(val_int1 * (bins+1) + val_int2) += ((binwidth1 - val_rem1) * (binwidth2 - val_rem2)) / binarea * wt(i);
    d_l(val_int1 * (bins+1) + val_int2) += ((binwidth1 - val_rem1) * (binwidth2 - val_rem2)) / binarea * Y(i) * wt(i);

    if(val_int2 < bins){
    c_l(val_int1 * (bins+1) + val_int2+1) += ((binwidth1 - val_rem1) * val_rem2) / binarea * wt(i);
    d_l(val_int1 * (bins+1) + val_int2+1) += ((binwidth1 - val_rem1) * val_rem2) / binarea * Y(i) * wt(i);
    }

    if (val_int1 < bins){
      c_l((val_int1+1) * (bins+1) + val_int2) += ((binwidth2 - val_rem2) * val_rem1) / binarea * wt(i);
      d_l((val_int1+1) * (bins+1) + val_int2) += ((binwidth2 - val_rem2) * val_rem1) / binarea * Y(i) * wt(i);
      if(val_int2 < bins){
      c_l((val_int1+1) * (bins+1) + val_int2+1) += (val_rem1 * val_rem2) / binarea * wt(i);
      d_l((val_int1+1) * (bins+1) + val_int2+1) += (val_rem1 * val_rem2) / binarea * Y(i) * wt(i);
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("weight") = c_l, 
                            Rcpp::Named("y") = d_l.array()/c_l.array(),
                            Rcpp::Named("x") = xgrid);
}

// [[Rcpp::export]]
Eigen::VectorXd llrt_cpp(const Eigen::MatrixXd &X, const Eigen::VectorXd &Y, const Eigen::MatrixXd& X_pred, const Eigen::VectorXd& wt, int method, int kcode,
                         double epsilon, const Eigen::VectorXd& h, int N_min){
    
    Eigen::VectorXd beta_0(X_pred.rows());
    kdtree tree(X, Y, wt, N_min);
    
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < X_pred.rows(); i++) {
      std::pair<Eigen::MatrixXd, Eigen::VectorXd> XtXXtY = tree.find_XtXXtY(X_pred.row(i), method, epsilon, h, kcode);
      beta_0(i) = calculate_beta(XtXXtY.first, XtXXtY.second);
    }
    return beta_0; 
} 


// [[Rcpp::export]]
Rcpp::List tgcv_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& Y, const Eigen::VectorXd& wt, int method, int kcode,
                        double epsilon, const Eigen::MatrixXd& bw, int N_min){

  Eigen::VectorXd bw_opt = Eigen::VectorXd::Zero(bw.cols());
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(X.cols()+1, X.cols()+1);
  Eigen::VectorXd SSEs(bw.rows());

  kdtree tree(X, Y, wt, N_min);
  double SSE_opt = -1;

  for (int i = 0; i< bw.rows(); i++) {            //calculating SSE for subsequent bandwidth
    Eigen::VectorXd h = bw.row(i);
    double SSE = 0;
    bool bwerror = false;
      
    #pragma omp parallel for reduction(+:SSE) schedule(static)
    for (int j = 0; j < X.rows(); j++) {
      double w = max_weight(kcode, h, wt(j)); 
      std::pair<Eigen::MatrixXd, Eigen::VectorXd> XtXXtY = tree.find_XtXXtY(X.row(j), method, epsilon, h, kcode);
      double beta_0 = calculate_beta(XtXXtY.first, XtXXtY.second);
      Eigen::MatrixXd XtX = XtXXtY.first;
      Eigen::MatrixXd XtXi = XtX.ldlt().solve(I);
      double w_point = w * XtXi(0,0);
//   Rcpp::Rcout << w_point << '\n';
      if (w_point == 1 || beta_0 == R_PosInf){
        bwerror = true; 
      }
      SSE += pow((beta_0 - Y(j))/(1-w_point), 2);
    }
    if (bwerror == true){
      SSEs(i) = R_PosInf;
    }
    else { 
      SSEs(i) = SSE;
    }
    
    if (SSE_opt == -1 && bwerror == false){
      SSE_opt = SSE;
      bw_opt = h;
    }
    else if (SSE <= SSE_opt && bwerror == false) {
      SSE_opt = SSE;
      bw_opt = h;
    }
  }
  return Rcpp::List::create(Rcpp::Named("bw_opt") = bw_opt, 
                            Rcpp::Named("SSE") = SSEs);
}


// [[Rcpp::export]]
Rcpp::List approx_gcv_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& Y, const Eigen::VectorXd& wt, int method, int kcode,
                    double epsilon, const Eigen::MatrixXd& bw, int N_min){
  
  int n = X.rows();
  Eigen::VectorXd bw_opt = Eigen::VectorXd::Zero(bw.cols());
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(X.cols()+1, X.cols()+1);
  Eigen::VectorXd SSE_clipped(n);
  Eigen::VectorXd SSEs(bw.rows());
  int cutoff = n * 0.95; 
  std::vector<std::pair<double, double>> w_SSE(n); 
  
  kdtree tree(X, Y, wt, N_min);
  double SSE_opt = -1;
  
  for (int i = 0; i< bw.rows(); i++) {            //calculating SSE for subsequent bandwidth
    Eigen::VectorXd h = bw.row(i);
    double SSE = 0;
    bool bwerror = false;
    
    #pragma omp parallel for reduction(+:SSE) schedule(static)
    for (int j = 0; j < n; j++) {
      double w = max_weight(kcode, h, wt(j)); 
      std::pair<Eigen::MatrixXd, Eigen::VectorXd> XtXXtY = tree.find_XtXXtY(X.row(j), method, epsilon, h, kcode);
      double beta_0 = calculate_beta(XtXXtY.first, XtXXtY.second);
      Eigen::MatrixXd XtX = XtXXtY.first;
      Eigen::MatrixXd XtXi = XtX.ldlt().solve(I);
      double w_point = w * XtXi(0,0);
      //   Rcpp::Rcout << w_point << '\n';
      if (w_point == 1 || beta_0 == R_PosInf){
        bwerror = true; 
      } 
      // double SSE = pow((beta_0 - Y(j))/(1-w_point), 2);
      // SSE_clipped(j) = SSE;
      w_SSE[j].first = w_point;
      w_SSE[j].second = beta_0-Y(j);
    }
    if (bwerror == true){
      SSEs(i) = R_PosInf;
    }
    else { 
//    std::nth_element(SSE_clipped.data(), SSE_clipped.data()+cutoff, SSE_clipped.data()+n);
      std::sort(w_SSE.begin(), w_SSE.end(), [](auto &left, auto &right) {
      return left.second < right.second;
      });
      for (int k = 0; k < cutoff; k++ ){
//    SSE_tot += SSE_clipped(k);
      SSE += pow(w_SSE[k].second/(1-w_SSE[k].first), 2);
        }
      SSEs(i) = SSE;
    }
    if (SSE_opt == -1 && bwerror == false){
      SSE_opt = SSE;
      bw_opt = h;
    }
    else if (SSE <= SSE_opt && bwerror == false) {
      SSE_opt = SSE;
      bw_opt = h;
    }
  }
  return Rcpp::List::create(Rcpp::Named("bw_opt") = bw_opt, 
                            Rcpp::Named("SSE") = SSEs);
}










