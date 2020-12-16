//#include <Eigen/Dense>
#include "RcppEigen.h"
#include "KDtreeiteration.h"
#include <utility>
#include <algorithm>
#include <iostream>   
#include <random>
#include <stack>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

Timer::Timer() : beg_(clock_::now()) {}

void Timer::reset() { beg_ = clock_::now(); }

double Timer::elapsed() const { 
  return std::chrono::duration_cast<second_>
  (clock_::now() - beg_).count(); }


kdtree::kdtree() = default; 
kdtree::~kdtree() = default; 

kdnode::kdnode() = default;  // constructor 
kdnode::kdnode(kdnode&& rhs) = default;  //move 
kdnode::~kdnode() = default;  // destructor  

template<typename T>            //overloading vector cout for printing purposes
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
  out << "{";
  size_t last = v.size() - 1;
  for(size_t i = 0; i < v.size(); ++i) {
    out << v[i];
    if (i != last) 
      out << ", ";
  }
  out << "}";
  return out;
}

std::unique_ptr<kdnode> kdtree::build_tree(all_point_t::iterator start, all_point_t::iterator end, 
                                           int split_d, double split_v, int N_min, size_t len,
                                           std::vector<double> dim_max_, std::vector<double> dim_min_){
  
  std::unique_ptr<kdnode> newnode = std::make_unique<kdnode>();
  
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
      Eigen::VectorXd XY = *i;    
      Eigen::VectorXd Y = XY.tail(1); 
      Eigen::VectorXd X = XY.head(XY.size()-1); 
      XtY_ += X*Y; 
      XtX_ += X*X.transpose();  
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


kdtree::kdtree(all_point_t& XY_arr, int N_min) {   // command to create kd tree 
  
  size_t len = XY_arr.size(); 
  
  std::vector<double> dim_max_;        
  std::vector<double> dim_min_; 
  
  for(size_t i=1; i<XY_arr[0].size()-1; i++){         
    dim_max_.push_back(XY_arr[0](i)); 
    dim_min_.push_back(XY_arr[0](i));
  }
  
  for(size_t i=1; i<XY_arr[0].size()-1; i++) { 
    for (size_t j=0; j<XY_arr.size(); j++){
      if (XY_arr[j](i) > dim_max_.at(i-1)){
        dim_max_.at(i-1) = XY_arr[j](i);
      }
      if (XY_arr[j][i] < dim_min_.at(i-1)){ 
        dim_min_.at(i-1) = XY_arr[j](i);
      }
    }
  } 
  
  root = build_tree(XY_arr.begin(), XY_arr.end(), 1, 1, N_min, len, dim_max_, dim_min_);   // using build_tree to build the entire tree 
}

all_point_t XYmat_to_XYarr(const Eigen::MatrixXd& XY_mat){    // conversion to vector form 
  std::vector<Eigen::VectorXd> XY_arr; 
  Eigen::MatrixXd XY_temp(XY_mat.rows(), XY_mat.cols()+1); 
  XY_temp.rightCols(XY_mat.cols()) = XY_mat; 
  XY_temp.col(0) = Eigen::VectorXd::Ones(XY_mat.rows()); 
  
  for (int i = 0; i <XY_temp.rows(); i++){
    XY_arr.push_back(XY_temp.row(i));
  }
  return XY_arr; 
}

inline double eval_kernel(int kcode, const double& z, double w)
{
  double tmp;
  if(abs(z) > 1) 
    return 0;
  else { 
    switch(kcode)
    {
    case 1: return 3*(1-z*z)/4*w; // Epanechnikov
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


inline Eigen::VectorXd get_w(const int kcode, const Eigen::VectorXd& x)
{
  Eigen::VectorXd temp(x.size());
  temp =  (1-x.array().pow(2)) * ((3./4));
  
  return temp;
}


std::pair<double,double> calculate_weight(int kcode, const Eigen::VectorXd& X_query, const std::vector<double>& dim_max, 
                                          const std::vector<double>& dim_min, const Eigen::VectorXd& h) {  // weight calculation of point
  if(dim_max.size()!= X_query.size()){
    //  std::cout << "error in dimensions"; 
    //    throw std::exception(); 
  }
  double w_max = 1;
  double w_min = 1; 
  double dis_max = 0; 
  double dis_min = 0;
  
  for(int i=0; i < dim_max.size(); i++) { 
    if (X_query(i) <= dim_max.at(i) && X_query(i) >= dim_min.at(i)) {
      dis_max = std::max(abs(dim_max.at(i) - X_query(i)), abs(dim_min.at(i) - X_query(i)));    // dis_max = which ever dim is further
      dis_min = 0;                                                                             
      w_min *= eval_kernel(kcode, dis_max/h(i)) / h(i);    // kern weight = multiplication of weight of each dim 
      w_max *= eval_kernel(kcode, dis_min/h(i)) / h(i);
      
    }
    else if (abs(X_query(i) - dim_max.at(i)) > abs(X_query(i) - dim_min.at(i))){
      dis_max = X_query(i) - dim_max.at(i);
      dis_min = X_query(i) - dim_min.at(i);
      w_min *= eval_kernel(kcode, dis_max/h(i)) / h(i);
      w_max *= eval_kernel(kcode, dis_min/h(i)) / h(i); 
    }
    else{
      dis_max = X_query(i) - dim_min.at(i);
      dis_min = X_query(i) - dim_max.at(i); 
      w_min *= eval_kernel(kcode, dis_max/h(i)) / h(i); 
      w_max *= eval_kernel(kcode, dis_min/h(i)) / h(i); 
    }
  }
  return std::make_pair(w_max, w_min); 
}

std::pair<Eigen::MatrixXd, Eigen::VectorXd> kdtree::get_XtXXtY(const Eigen::VectorXd& X_query,  // to obtain the exact XtXXtY from tree 
                                                               const std::vector<double>& dim_max, 
                                                               const std::vector<double>& dim_min, 
                                                               std::unique_ptr<kdnode>& root, 
                                                               const Eigen::VectorXd& h, 
                                                               int kcode){
  
  std::pair<double,double> w_maxmin;
  //std::stack<std::unique_ptr<kdnode>> storage; 
  std::stack <kdnode *> storage; 
  //std::unique_ptr<kdnode> curr = root; 
  auto curr = root.get(); 
  Eigen::MatrixXd XtX = Eigen::MatrixXd::Zero(curr->XtX.rows() , curr->XtX.cols()); 
  Eigen::VectorXd XtY = Eigen::MatrixXd::Zero(curr->XtY.rows() , curr->XtY.cols());
  w_maxmin = calculate_weight(kcode, X_query, dim_max, dim_min,h);
  double w_max = w_maxmin.first; 
  double w_min = w_maxmin.second;
  
  while (w_max != w_min || storage.empty() == false){
    while (w_max != w_min ){   
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
  return std::make_pair(XtX,XtY); 
}

std::pair<Eigen::MatrixXd, Eigen::VectorXd> kdtree::getapprox_XtXXtY(const Eigen::VectorXd& X_query,    // to obtain the approx XtXXtY from tree 
                                                                     const std::vector<double>& dim_max,
                                                                     const std::vector<double>& dim_min, 
                                                                     std::unique_ptr<kdnode>& root,
                                                                     double epsilon, 
                                                                     const Eigen::VectorXd& h,
                                                                     int kcode){ 
  
  std::pair<double,double> w_maxmin;
  //std::stack<std::unique_ptr<kdnode>> storage; 
  //std::unique_ptr<kdnode> curr = root; 
  std::stack<kdnode*> storage; 
  auto curr =root.get(); 
  Eigen::MatrixXd XtX = Eigen::MatrixXd::Zero(curr->XtX.rows() , curr->XtX.cols()); 
  Eigen::VectorXd XtY = Eigen::MatrixXd::Zero(curr->XtY.rows() , curr->XtY.cols());
  
  w_maxmin = calculate_weight(kcode, X_query, dim_max, dim_min, h);
  double w_max = w_maxmin.first; 
  double w_min = w_maxmin.second;
  double weight_sf = 0; 
  
  while (( curr->left_child.get() != nullptr && w_max-w_min > 2*epsilon*(weight_sf + curr->n_below*w_min)) || storage.empty() == false  ){
    while (w_max - w_min > 2*epsilon*(weight_sf + curr->n_below*w_min) && curr->left_child.get() != nullptr ){   // if condition fufilled        
      storage.push(curr);
      curr = curr->left_child.get();  
      w_maxmin = calculate_weight(kcode, X_query, curr->dim_max, curr->dim_min, h);  // calculate max and min weight 
      w_max = w_maxmin.first;                      
      w_min = w_maxmin.second; 
      if(w_max - w_min <= 2 * epsilon * (weight_sf + curr->n_below*w_min)){ 
        double weight = 0.5*(w_max + w_min);
        weight_sf += weight;  
        XtX += weight * curr->XtX;
        XtY += weight * curr->XtY; 
      }           
    }
    
    curr = storage.top(); 
    storage.pop(); 
    curr = curr->right_child.get(); 
    w_maxmin = calculate_weight(kcode, X_query, curr->dim_max, curr->dim_min, h);  // calculate max and min weight 
    w_max = w_maxmin.first;   
    w_min = w_maxmin.second; 
    if(w_max - w_min <= 2 * epsilon * (weight_sf + curr->n_below*w_min)){ 
      double weight = 0.5*(w_max + w_min);
      weight_sf += weight;  
      XtX += weight * curr->XtX;
      XtY += weight * curr->XtY; 
    }          
  }
  return std::make_pair(XtX,XtY); 
}

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


Eigen::VectorXd form_ll_XtY(const Eigen::VectorXd& XtY, const Eigen::VectorXd& X_query){ //steps to form the local linear XtY(with query pts)
  Eigen::VectorXd extra_XtY = Eigen::VectorXd::Zero((XtY.size()));
  Eigen::VectorXd ll_XtY(XtY.size());
  extra_XtY.tail(XtY.size()-1) = X_query * XtY.head(1); 
  ll_XtY = XtY - extra_XtY; 
  return ll_XtY;
}

Eigen::VectorXd calculate_beta(int kcode, const Eigen::MatrixXd &XtX, const Eigen::VectorXd &XtY) {
  Eigen::LDLT<Eigen::MatrixXd> LDLT_XtX(XtX);
  return(LDLT_XtX.solve(XtY));;  //first term gives the predicted value
}

std::pair<Eigen::VectorXd, double> calculate_beta_XtXinv(int kcode, const Eigen::MatrixXd &XtX, 
                                                         const Eigen::MatrixXd &XtY) { 
  
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(XtX.rows(),XtX.cols()); 
  Eigen::MatrixXd XtX_inv = XtX.ldlt().solve(I);          
  
  return std::make_pair(XtX.ldlt().solve(XtY), XtX_inv(0,0)); //first term gives the predicted value, second term gives the term used for loocv 
}

double max_weight(int kcode, const Eigen::VectorXd&h){ //return max_weight possible of a point 
  double w_max = 1; 
  for(size_t i =0; i < h.size(); i++){
    w_max *= eval_kernel(kcode, 0)/h(i);
  }
  return w_max;
}

//' Multiplies two doubles
//'
//' @param X X 
//' @param Y Y 
//' @param X_pred X_pred
//' @param kcode kcode
//' @param h h 
//' @export
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
       double w = eval_kernel(kcode, dif/h, wt(j));
       s0 += w;
       s1 += w * dif;
       s2 += w * pow(dif,2);
       t0 += w * Y(j+lb);
       t1 += w * dif * Y(j+lb);
    }
     double beta = (s2*t0 - s1*t1)/(s0*s2 - pow(s1,2)); 
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
        double w = eval_kernel(kcode,x1/h(0), wt(j)) * eval_kernel(kcode, x2/h(1), wt(j));
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
    Eigen::VectorXd beta = calculate_beta(kcode, s, t);
    beta_0(i) = beta(0);
  }
  return beta_0; 
}


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
      double w = 1; 
      for (int k = 0; k < X.cols(); k++) {
        w *= eval_kernel(kcode, Xdif(k)/h(k), wt(j)); 
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
    Eigen::VectorXd beta = calculate_beta(kcode, s, t);
    beta_0(i) = beta(0);
  }
  return beta_0; 
}

// [[Rcpp::export]]
Rcpp::List bin1d_cpp(const Eigen::VectorXd &X, const Eigen::VectorXd&Y, int bins){
  
  double grid_start = X.minCoeff(); 
  double grid_end = X.maxCoeff(); 
  double binwidth = (grid_end - grid_start) / (bins);
  
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
    c_l(val_int) += (1 - val_rem); 
    d_l(val_int) += (1 - val_rem) * Y(i);
    
    if (val_int < bins) {
      c_l(val_int + 1) += val_rem ; 
      d_l(val_int + 1) += val_rem * Y(i); 
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("xweight") = c_l, 
                            Rcpp::Named("y") = d_l.array()/c_l.array(),
                            Rcpp::Named("x") = xgrid);
}
  
//   std::vector<double> kern_w;
//   
//   int L_i = int(round(h/binwidth)); // bins-1); 
//   for(int i = kern_w.size(); i <= L_i; i++){ 
//     double w = eval_kernel(kcode, (i*binwidth)/h)/ h;
//     kern_w.push_back(w);
//   }
//   
//   Eigen::VectorXd results(X_pred.size()); 
//   for(int i = 0; i < X_pred.size(); i++ ){
//     double s_0 = 0; 
//     double s_1 = 0; 
//     double s_2 = 0; 
//     double t_0 = 0; 
//     double t_1 = 0;  
//     
//     int Xpred_bin = round((X_pred(i) - grid_start) / binwidth); 
//     int min = std::max(Xpred_bin - L_i, 0);
//     int max = std::min(Xpred_bin + L_i, bins-1);
//     
//     for(int j = min ; j <= max; j++) {
//       int bin_dif = abs(Xpred_bin - j); 
//       
//       s_0 += kern_w[bin_dif] * c_l(j);
//       s_1 += kern_w[bin_dif] * (bin_dif * binwidth) * c_l(j);
//       s_2 += kern_w[bin_dif] * pow((bin_dif * binwidth), 2) * c_l(j); 
//       t_0 += kern_w[bin_dif] * d_l(j); 
//       t_1 += kern_w[bin_dif] * (bin_dif * binwidth) * d_l(j); 
//       
//     }
//     results(i) = ((s_2 * t_0) - (s_1*t_1))/((s_2 * s_0)-pow(s_1, 2));
//   }  
//   return results; 
// }


// Eigen::VectorXd loclin(const Eigen::MatrixXd& X, const Eigen::VectorXd& Y, int method, int kcode, 
//                        double epsilon, const Eigen::VectorXd& h, int N_min){    //method 1 for exact, 2 for approximate, bandwidth is 1d, will add vector later, epsilon is for epsilon for
//   
//   Eigen::MatrixXd XY_mat(X.rows(), X.cols() + 1);
//   XY_mat << X, Y; 
//   
//   all_point_t XY_arr = XYmat_to_XYarr(XY_mat);   //preprocessing work to convert data into structures that we want 
//   
//   Eigen::VectorXd beta_0(XY_arr.size());
//   std::pair<Eigen::MatrixXd, Eigen::VectorXd> XtXXtY;
//   Eigen::VectorXd beta;
//   
//   kdtree tree(XY_arr, N_min);     //building tree 
//   
//   #pragma omp parallel for schedule(static) private(XtXXtY, beta)
//   for (size_t i = 0; i< XY_arr.size(); i++) { 
//     XtXXtY = tree.find_XtXXtY(X.row(i), method, epsilon, h, kcode); 
//     beta = calculate_beta(kcode, XtXXtY.first, XtXXtY.second);
//     beta_0(i) = beta(0); 
//   }
//   return beta_0; 
// }

// [[Rcpp::export]]
Eigen::VectorXd predict(const Eigen::MatrixXd& X, const Eigen::VectorXd& Y, const Eigen::MatrixXd& X_pred, int method, int kcode, 
                        double epsilon, const Eigen::VectorXd& h,  int N_min){    //method 1 for exact, 2 for approximate, bandwidth is 1d, will add vector later, epsilon is for epsilon for
  
  Eigen::MatrixXd XY_mat(X.rows(), X.cols() + 1);
  XY_mat << X, Y; 
  
  all_point_t XY_arr = XYmat_to_XYarr(XY_mat) ;     //preprocessing work to convert data into structures that we want 
  Eigen::VectorXd beta_0(X_pred.rows());
  
  kdtree tree(XY_arr, N_min); //building tree
  
  #pragma omp parallel for schedule(static) 
  for (size_t i = 0; i< X_pred.rows(); i++) { 
    std::pair<Eigen::MatrixXd, Eigen::VectorXd> XtXXtY;
    XtXXtY = tree.find_XtXXtY(X_pred.row(i), method, epsilon, h, kcode); 
    Eigen::VectorXd beta = calculate_beta(kcode, XtXXtY.first, XtXXtY.second);
    beta_0(i) = beta(0); 
  }
  return beta_0; 
}

// [[Rcpp::export]]
Eigen::VectorXd bw_loocv(const Eigen::MatrixXd& X, const Eigen::VectorXd &Y, int method, int kcode, double epsilon,
                         const Eigen::MatrixXd& bw, int N_min){
  
  Eigen::MatrixXd XY_mat(X.rows(), X.cols() + 1);
  XY_mat << X, Y;  
  
  Eigen::VectorXd SSE_arr(bw.rows());
  Eigen::VectorXd bw_opt;
  
  all_point_t XY_arr = XYmat_to_XYarr(XY_mat);
  Eigen::VectorXd beta_0(XY_arr.size());
  
  kdtree tree(XY_arr, N_min);
  
  double SSE_opt = -1;
  double SSE;
  double w_point;
  
  Eigen::VectorXd h = bw.row(0);
  
  for (int i = 0; i< bw.rows(); i++) {            //calculating SSE for subsequent bandwidth
    h = bw.row(i);
    SSE = 0;
    
    #pragma omp parallel for reduction(+: SSE)
    for (int j = 0; j< XY_arr.size(); j++) {
      std::pair<Eigen::MatrixXd, Eigen::VectorXd> XtXXtY;
      XtXXtY = tree.find_XtXXtY(X.row(j), method, epsilon, h, kcode);
      std::pair<Eigen::VectorXd, double> beta_XtXinv = calculate_beta_XtXinv(kcode, XtXXtY.first, XtXXtY.second);
      Eigen::VectorXd beta = beta_XtXinv.first;
      beta_0(j) = beta(0);
      w_point = max_weight(kcode, h) * beta_XtXinv.second;
      SSE += pow((beta_0(j) - Y(j)) / (1-w_point),2);
    }
    SSE_arr(i) = SSE;
    if (SSE_opt == -1 && SSE >= 0){ 
      SSE_opt = SSE;
      bw_opt = h;
    }
    else if (SSE <= SSE_opt) {
      SSE_opt = SSE;
      bw_opt = h;
    }
  }
  return bw_opt;
}
