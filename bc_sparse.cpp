// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
SEXP sparseBC( SEXP XX_ ) {
  typedef Eigen::MatrixXd Mat;  
  typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMat;
  typedef Eigen::SparseVector<double, Eigen::RowMajor> SpVec;
  typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
  MapMatd X(Rcpp::as<MapMatd>(XX_));
  SpMat Xsparse = X.sparseView();
  Mat res = Eigen::MatrixXd::Zero(X.rows(),X.rows());
  SpVec xi(X.cols()), xj(X.cols());

  double Cij, Si_plus_Sj, min;
  
  for(int i=0; i<X.rows(); i++){
    for(int j=0; j<X.rows(); j++){
      if(i != j){
        xi = Xsparse.row(i);
        xj = Xsparse.row(j);
        Cij = xi.cwiseMin(xj).sum();
        Si_plus_Sj = xi.sum() + xj.sum();
        res(i, j) = (Si_plus_Sj - 2.0 * Cij ) / Si_plus_Sj; 
      } else {
        res(i, j) = 0.0;
      }
        
      
    }
  }
return(wrap(res));

}

/*** R
my_bray_curtis <- function(df){
  dis <- sparseBC(as.matrix(df))
  colnames(dis) <- row.names(df)
  row.names(dis) <- row.names(df)
  return(dis)
}
  */



