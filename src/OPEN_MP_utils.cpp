#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]

//' Get number of processors
//' @return number of processors 
//' @export
// [[Rcpp::export]]
int get_num_procs() {
  #ifdef _OPENMP
    return (omp_get_num_procs());
  #else
    return (1);
  #endif
}

//' Set number of threads to use 
//' @param threads number of threads to use 
//' @return number of threads set 
//' @export
// [[Rcpp::export]]
int set_num_threads(int threads) {
#ifdef _OPENMP
  omp_set_num_threads(threads);
  return (omp_get_max_threads());
#else
  return (1L);
#endif
}