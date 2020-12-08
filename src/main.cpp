// #include <iostream>
// //#include "KDtree.h"
// #include "KDtreeiteration.h"
// #include <utility>
// #include <iostream>
// //#include "KDtreetest.h"
// 
// // Eigen::VectorXd loclin(const Eigen::MatrixXd& XY_mat, int method, int kcode, 
//                          // double epsilon, const Eigen::VectorXd& h, int N_min);
// int main(){ 
//    Eigen::MatrixXd test1 = Eigen::MatrixXd::Random(1000,2);
//    Eigen::MatrixXd test2 = Eigen::VectorXd::Random(1000);
//    Eigen::VectorXd test3 = Eigen::VectorXd::Random(200); 
//    Eigen::VectorXd h(2); 
//    h << 0.2, 0.1;
//    // h<< 0.2;
//    // Timer tmr; 
//    // predict1d(test1, test2, 1, 0.2); 
//    // double t = tmr.elapsed() ; 
//    // std::cout << t << '\n'; 
//    // tmr.reset();
//    // predict(test1, test2, 1, 1,0, h, 1);
//    // t = tmr.elapsed();
//    // std::cout << t;
//    //h << 0.2;
//    //std::cout << bin1d(test1, test2, 1, 0.2, 1000) - predict1d(test1, test2, 1, 0.2);
//    //std::cout << bin1d(test1, test2, 1, 0.2 , 800);
//    //Eigen::VectorXd z  = test1.col(0); 
//    //std::cout <<"sum" << z.sum() <<'\n';
//    //predict2dd(test1, test2, test1, 1, h);
//    loclin(test1, test2, 1 , 1, 0, h, 1);
//    
// }  