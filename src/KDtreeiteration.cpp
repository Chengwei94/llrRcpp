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

std::shared_ptr<kdnode> kdtree::build_tree(all_point_t::iterator start, all_point_t::iterator end, 
                                           int split_d, double split_v, int N_min, size_t len,
                                           const std::vector<double>& dim_max_, const std::vector<double>& dim_min_){
    
    std::shared_ptr<kdnode> newnode = std::make_shared<kdnode>();
    
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
            if (XY_arr[j](i) > dim_max.at(i-1)){
                dim_max_.at(i-1) = XY_arr[j](i);
            }
            if (XY_arr[j][i] < dim_min.at(i-1)){ 
                dim_min_.at(i-1) = XY_arr[j](i);
            }
        }
    } 
    
    root = build_tree(XY_arr.begin(), XY_arr.end(), 1, 1, N_min, len, dim_max, dim_min);   // using build_tree to build the entire tree 
}

all_point_t XYmat_to_Xarr(const Eigen::MatrixXd& XY_mat){      
    Eigen::MatrixXd XY_trans = XY_mat.transpose(); 
    all_point_t X_arr;
    for (size_t i=0; i<XY_trans.cols(); i++) {
        Eigen::VectorXd x; 
        x = XY_trans.block(0, i, XY_trans.rows()-1, 1); 
        X_arr.push_back(x);
    }
    return X_arr; 
}

all_point_t Xmat_to_Xarr(const Eigen::MatrixXd& X_mat){     //  conversion to X_query form 
    std::vector<Eigen::VectorXd> X_arr; 
    Eigen::MatrixXd X_trans = X_mat.transpose(); 
    for (size_t i=0; i<X_trans.cols(); i++) {
        Eigen::VectorXd x; 
        x = X_trans.block(0, i, X_trans.rows(), 1); 
        X_arr.push_back(x);
    }
    return X_arr; 
}

all_point_t XYmat_to_XYarr(const Eigen::MatrixXd& XY_mat){    // conversion to vector form 
    std::vector<Eigen::VectorXd> XY_arr; 
    Eigen::MatrixXd XY_temp(XY_mat.rows(), XY_mat.cols()+1); 
    XY_temp.topRightCorner(XY_mat.rows(), XY_mat.cols()) = XY_mat; 
    XY_temp.col(0) = Eigen::VectorXd::Ones(XY_mat.rows()); 
    
    for (int i = 0; i <XY_temp.rows(); i++){
        XY_arr.push_back(XY_temp.row(i));
    }
    return XY_arr; 
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
                                                               std::vector<double> dim_max, 
                                                               std::vector<double> dim_min, 
                                                               std::shared_ptr<kdnode>& root, 
                                                               const Eigen::VectorXd& h, 
                                                               int kcode){
    
    std::pair<double,double> w_maxmin;
    std::stack<std::shared_ptr<kdnode>> storage; 
    std::shared_ptr<kdnode> curr = root; 
    Eigen::MatrixXd XtX = Eigen::MatrixXd::Zero(curr->XtX.rows() , curr->XtX.cols()); 
    Eigen::VectorXd XtY = Eigen::MatrixXd::Zero(curr->XtY.rows() , curr->XtY.cols());
    w_maxmin = calculate_weight(kcode, X_query, dim_max, dim_min,h);
    double w_max = w_maxmin.first; 
    double w_min = w_maxmin.second;
    
    while (w_max != w_min || storage.empty() == false){
        while (w_max != w_min ){   
            storage.push(curr);
            curr = curr->left_child;  
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
        curr = curr->right_child; 
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
                                                                     std::vector<double> dim_max,
                                                                     std::vector<double> dim_min, 
                                                                     std::shared_ptr<kdnode>& root,
                                                                     double epsilon, 
                                                                     const Eigen::VectorXd& h,
                                                                     int kcode){ 
    
    std::pair<double,double> w_maxmin;
    std::stack<std::shared_ptr<kdnode>> storage; 
    std::shared_ptr<kdnode> curr = root; 
    Eigen::MatrixXd XtX = Eigen::MatrixXd::Zero(curr->XtX.rows() , curr->XtX.cols()); 
    Eigen::VectorXd XtY = Eigen::MatrixXd::Zero(curr->XtY.rows() , curr->XtY.cols());
    
    w_maxmin = calculate_weight(kcode, X_query, dim_max, dim_min, h);
    double w_max = w_maxmin.first; 
    double w_min = w_maxmin.second;
    double weight_sf = 0; 
    
    while (( curr->left_child != nullptr && w_max-w_min > 2*epsilon*(weight_sf + curr->n_below*w_min)) || storage.empty() == false  ){
        while (w_max - w_min > 2*epsilon*(weight_sf + curr->n_below*w_min) && curr->left_child != nullptr ){   // if condition fufilled        
            storage.push(curr);
            curr = curr->left_child;  
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
        curr = curr->right_child; 
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

// [[Rcpp::export]]
Eigen::VectorXd predict1d(const Eigen::MatrixXd& XY_mat, Eigen::VectorXd& X_pred, int kcode, double h){ 
    //all_point_t XY_arr = XYmat_to_XYarr(XY_mat);   //preprocessing work to convert data into structures that we want 
    //std::sort(XY_arr.begin(), XY_arr.end(), [](const Eigen::VectorXd& a, const Eigen::VectorXd & b) {
    //    return a(1) < b(1); 
    //});  
    //std::sort(X_pred.data(), X_pred.data() + X_pred.size());
    
    //int ssize = XY_arr.size(); // ssize = sample size
    // now XY_arr is sorted, so is X_pred: 
    //std::vector<Eigen::MatrixXd> XtX(ssize, Eigen::MatrixXd(2,2)); 
    //std::vector<Eigen::VectorXd> XtY(ssize, Eigen::VectorXd(2));
    //std::vector<int> checker(XY_arr.size(), -1);  // creating a checker of -1; 
    Eigen::VectorXd beta_0(X_pred.size());
    //std::cout << XY_arr << '\n'; 
    //std::cout << X_pred << '\n'; 
    
    int lb = 0; 
    int ub = 0; 
    Eigen::VectorXd X(2);
    X(0) = 1;
    
    for(int i= 0; i < X_pred.size(); i++) { 
        while (lb < XY_mat.rows() && XY_mat(lb, 0)< (X_pred(i) - h)){
            lb++; 
        } 
        while (ub < XY_mat.rows() && XY_mat(ub ,0) < (X_pred(i) + h)){ 
            ub++;
        }
        Eigen::MatrixXd XtWX(2,2); 
        XtWX.setZero();
        Eigen::VectorXd XtWY(2);
        XtWY.setZero(); 
        for(int j = lb; j < ub; j++){ 
            //std::cout << XtWY << '\n';
            X(1) = XY_mat(j, 0) - X_pred(i); 
            //Eigen::MatrixXd XtX = X*X.transpose(); 
            double Y = XY_mat(j, 1);
           // Eigen::MatrixXd XtY = X*Y;
            
            double w = eval_kernel(kcode, abs(X_pred(i) - XY_mat(j, 0))/h)/h;
            //std::cout << "w[ij]"  << w << i << j << '\n' ;
            //XtWX += w*XtX; 
            //XtWY += w*XtY;
        }
        //std::cout << Xt << "\n"; 
        //Eigen::VectorXd beta = calculate_beta(kcode, XtWX, XtWY);
        //beta_0(i) = beta(0);
    }
    return beta_0;                 
}

// [[Rcpp::export]]
Eigen::VectorXd multi (Eigen::VectorXd &Z, Eigen::VectorXd &Y, Eigen::VectorXd& X_pred, int kcode, double h){ 
    Eigen::MatrixXd X(Z.size(), 2); 
    Eigen::VectorXd beta_0(X_pred.size());
    Eigen::VectorXd W(Z.size()); 
    for (int i = 0; i < X_pred.size(); i++){
        X.setOnes(); 
        for (int j = 0; j < Z.size() ; j++ ){ 
            X(j,1) = Z(j) - X_pred(i); 
        }
        // X.col(1) = Z.array() - X_pred(i); 
        W = get_w(kcode, X.col(1)/h);
        Eigen::VectorXd beta = calculate_beta(kcode, X.transpose()*W.asDiagonal()*X, X.transpose()*W.asDiagonal()*Y);
        beta_0(i) = beta(0); 
    }
    return beta_0; 
}

// [[Rcpp::export]]
Eigen::VectorXd predict1dd (const Eigen::VectorXd &X, const Eigen::VectorXd &Y, const Eigen::VectorXd& X_pred, int kcode, double h){
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
        Eigen::MatrixXd X1(m, 2);
        X1.setOnes(); 
        Eigen::VectorXd Y1(m);
        Eigen::VectorXd W(m); 
        for (int j = 0; j < m; j ++){
            X1(j, 1) = X(j+lb) - X_pred(i); 
            W(j) = eval_kernel(kcode,X1(j,1));
            Y1(j) = Y(j+lb); 
        }
        // Eigen::VectorXd W = get_w(kcode, X1.col(1)/h); 
    
        Eigen::VectorXd beta = calculate_beta(kcode, X1.transpose()*W.asDiagonal()*X1 , X1.transpose() * W.asDiagonal() *Y1);
        beta_0(i) = beta(0);
    }
    return beta_0; 
}

// [[Rcpp::export]]
Eigen::VectorXd predict2dd (const Eigen::MatrixXd &X, const Eigen::VectorXd &Y, const Eigen::MatrixXd& X_pred, int kcode, Eigen::VectorXd h){
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
        
        std::vector<int> idx;
        for (int j = lb; j < ub; j++) { 
            if (abs(X(j,1) - X_pred(i,1)) <= h(1))  
            idx.push_back(j);  
        }
        Eigen::MatrixXd X1(idx.size(), 3); 
        X1.setOnes(); 
        Eigen::VectorXd Y1(idx.size());
        Eigen::VectorXd W(idx.size());
        
        if(idx.size() <= 5){ 
          Rcpp::stop("Choose a larger bandwidth");
        }
        
        for (int j = 0; j < idx.size(); j ++){
            X1(j,1) = X(idx[j],0) - X_pred(i,0);
            X1(j,2) = X(idx[j],1) - X_pred(i,1);
            W(j)  = eval_kernel(kcode, X1(j,1)/h(0)) * eval_kernel (kcode, X1(j,2)/h(1));  
            Y1(j) = Y(idx[j]); 
        }

        // Eigen::VectorXd W = get_w(kcode, X1.col(1)/h); 
        Eigen::VectorXd beta = calculate_beta(kcode, X1.transpose()*W.asDiagonal()*X1 , X1.transpose() * W.asDiagonal() *Y1);
        beta_0(i) = beta(0);
    }
    return beta_0; 
}


// [[Rcpp::export]]
Eigen::VectorXd bin1d(const Eigen::VectorXd &X, const Eigen::VectorXd &Y, const Eigen::VectorXd& X_pred, int kcode, double h, int bins){
    
    double grid_start = X.minCoeff(); 
    double grid_end = X.maxCoeff(); 
    double binwidth = (grid_end - grid_start) / (bins-1);
   
    //if (X.size() < bins || h < binwidth){ 
      //  Rcpp::Rcout << "bandwidth too small or datapoints > bins, using non-bin local linear instead."; 
     //   return predict1dd(X, Y, X_pred, kcode, h);  
    //}

    // Linear binning;
    Eigen::VectorXd c_l = Eigen::VectorXd::Zero(bins);  
    Eigen::VectorXd d_l = Eigen::VectorXd::Zero(bins); 
    // add weights; 
    for (int i=0; i < X.size() ; i++) {
        double val = (X(i) - grid_start)/binwidth;
        int val_int = floor(val); 
        double val_rem = val - val_int; 
        c_l(val_int) += (1 - val_rem); 
        d_l(val_int) += (1 - val_rem) * Y(i);
        
        if (val_int < bins-1 ) {
            c_l(val_int + 1) += val_rem ; 
            d_l(val_int + 1) += val_rem * Y(i); 
        }
    }

    std::vector<double> kern_w;
    
    int L_i = int(round(h/binwidth)); // bins-1); 
    int L_i_neg = 0 - L_i;  
    for(int i = kern_w.size(); i <= L_i; i++){ 
        double w = eval_kernel(kcode, (i*binwidth)/h)/ h;
        kern_w.push_back(w);
    }
    
    Eigen::VectorXd results(X_pred.size()); 
    for(int i = 0; i < X_pred.size(); i++ ){
        double s_0 = 0; 
        double s_1 = 0; 
        double s_2 = 0; 
        double t_0 = 0; 
        double t_1 = 0;  
        
        int Xpred_bin = round((X_pred(i) - grid_start) / binwidth); 
        int min = std::max(Xpred_bin - L_i, 0);
        int max = std::min(Xpred_bin + L_i, bins-1);
        
        for(int j = min ; j <= max; j++) {
            int bin_dif = abs(Xpred_bin - j); 

            s_0 += kern_w[bin_dif] * c_l(j);
            s_1 += kern_w[bin_dif] * (bin_dif * binwidth) * c_l(j);
            s_2 += kern_w[bin_dif] * pow((bin_dif * binwidth), 2) * c_l(j); 
            t_0 += kern_w[bin_dif] * d_l(j); 
            t_1 += kern_w[bin_dif] * (bin_dif * binwidth) * d_l(j); 
            
        }
        results(i) = ((s_2 * t_0) - (s_1*t_1))/((s_2 * s_0)-pow(s_1, 2));
    }  
    return results; 
}

// [[Rcpp::export]]
Eigen::VectorXd loclin(const Eigen::MatrixXd& X, const Eigen::VectorXd& Y, int method, int kcode, 
                       double epsilon, const Eigen::VectorXd& h, int N_min){    //method 1 for exact, 2 for approximate, bandwidth is 1d, will add vector later, epsilon is for epsilon for
    
    Eigen::MatrixXd XY_mat(X.rows(), X.cols() + 1);
    XY_mat << X, Y; 

    all_point_t XY_arr = XYmat_to_XYarr(XY_mat);   //preprocessing work to convert data into structures that we want 

    Eigen::VectorXd beta_0(XY_arr.size());
    std::pair<Eigen::MatrixXd, Eigen::VectorXd> XtXXtY;
    Eigen::VectorXd beta;
   
    kdtree tree(XY_arr, N_min);     //building tree 
    
   // #pragma omp parallel for schedule(static) private(XtXXtY, beta)
    for (size_t i = 0; i< XY_arr.size(); i++) { 
        XtXXtY = tree.find_XtXXtY(X.row(i), method, epsilon, h, kcode); 
        beta = calculate_beta(kcode, XtXXtY.first, XtXXtY.second);
        beta_0(i) = beta(0); 
    }
    return beta_0; 
}

// [[Rcpp::export]]
Eigen::VectorXd predict(const Eigen::MatrixXd& X, const Eigen::VectorXd& Y, const Eigen::MatrixXd& X_pred, int method, int kcode, 
                        double epsilon, const Eigen::VectorXd& h,  int N_min){    //method 1 for exact, 2 for approximate, bandwidth is 1d, will add vector later, epsilon is for epsilon for
    
    Eigen::MatrixXd XY_mat(X.rows(), X.cols() + 1);
    XY_mat << X, Y; 
    
    all_point_t XY_arr = XYmat_to_XYarr(XY_mat) ;     //preprocessing work to convert data into structures that we want 
    Eigen::VectorXd beta_0(X_pred.rows());
    
    kdtree tree(XY_arr, N_min); //building tree
    
    // #pragma omp parallel for schedule(static) 
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
        
        // #pragma omp parallel for reduction(+: SSE)
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





    
    

