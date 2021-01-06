#pragma once
#include "RcppEigen.h"
#include <memory> 
#include <vector>
#include <chrono>
#include <iostream>


all_point_t XYmat_to_XYarr(const Eigen::MatrixXd& XY_mat);  
double eval_kernel(int kcode, const double& z);
double eval_kernel(int kcode, const double& z, double w);
std::pair<double, double> calculate_weight(int kcode, const Eigen::VectorXd& X_query, const std::vector<double>& dim_max, const std::vector<double>& dim_min , const Eigen::VectorXd& h); 
Eigen::MatrixXd form_ll_XtX(const Eigen::MatrixXd& XtX, const Eigen::VectorXd& X_query ); 
Eigen::VectorXd form_ll_XtY(const Eigen::VectorXd& XtY , const Eigen::VectorXd& X_query);
Eigen::VectorXd calculate_beta(const Eigen::MatrixXd& XtX, const Eigen::VectorXd& XtY);
std::pair<Eigen::VectorXd, double> calculate_beta_XtXinv(int kcode, const Eigen::MatrixXd &XtX, const Eigen::MatrixXd &XtY);
double max_weight(int kcode, const Eigen::VectorXd&h);

// R function
// Eigen::VectorXd loclin(const Eigen::MatrixXd& XY_mat, int method, int kcode, 
//                        double epsilon, const Eigen::VectorXd& h, int N_min);
Rcpp::List bin1d_cpp(const Eigen::VectorXd& X, const Eigen::VectorXd&Y, int bins, const Eigen::VectorXd& wt);
Eigen::VectorXd llr1d_cpp(const Eigen::VectorXd& X, const Eigen::VectorXd& Y, const Eigen::VectorXd& X_pred, int kcode, double h, Eigen::VectorXd wt);
Eigen::VectorXd llr2d_cpp (const Eigen::MatrixXd &X, const Eigen::VectorXd &Y, const Eigen::MatrixXd& X_pred, int kcode, Eigen::VectorXd h, Eigen::VectorXd wt);
Eigen::VectorXd llr_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& Y, const Eigen::MatrixXd& X_pred,
                        int kcode, Eigen::VectorXd h, const Eigen::VectorXd &wt);
Rcpp::List bin2d_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd&Y, int bins, const Eigen::VectorXd& wt);// bins at each dir
Eigen::VectorXd llrt_cpp(const kdtree& tree, const Eigen::MatrixXd& xpred, int method, int kcode, 
                         double epsilon, const Eigen::VectorXd& h, int N_min);