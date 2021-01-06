#ifndef kdtree_h
#define kdtree_h

#include "RcppEigen.h"
#include <memory> 
#include <vector>
#include <chrono>
#include <iostream>

using all_point_t = std::vector<Eigen::VectorXd>; 

class kdtree{
public:    
    // constructor
    kdtree(const Eigen::MatrixXd& X, const Eigen::VectorXd& Y, const Eigen::VectorXd& w, int N_min);
    
    class kdnode{ 
    public: 
        int n_below; 
        int split_d; 
        double split_v;
        Eigen::MatrixXd XtX; 
        Eigen::VectorXd XtY; 
        std::vector<double> dim_max; 
        std::vector<double> dim_min;
        std::unique_ptr<kdnode> right_child; 
        std::unique_ptr<kdnode> left_child;
    };
    std::unique_ptr<kdtree::kdnode> root; 
    std::unique_ptr<kdtree::kdnode> build_tree(all_point_t::iterator, all_point_t::iterator, int, double, int, size_t, std::vector<double> , std::vector<double> );    std::pair<Eigen::MatrixXd, Eigen::VectorXd> get_XtXXtY(const Eigen::VectorXd& X_query, const std::vector<double>& dim_max, const std::vector<double>& dim_min, std::unique_ptr<kdtree::kdnode>& root, const Eigen::VectorXd& h, int kcode);
    std::pair<Eigen::MatrixXd, Eigen::VectorXd> getapprox_XtXXtY(const Eigen::VectorXd& X_query, const std::vector<double>& dim_max, const std::vector<double>& dim_min, std::unique_ptr<kdtree::kdnode>& root, double epsilon, const Eigen::VectorXd& h, int kcode);
    std::pair<Eigen::MatrixXd, Eigen::VectorXd> find_XtXXtY(const Eigen::VectorXd& X_query, int method, double epsilon, const Eigen::VectorXd& h, int kcode);
};

#endif

