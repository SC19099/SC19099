#include <Rcpp.h>
using namespace Rcpp;

double vacc3a(double age, bool female, bool ily){
  double p = 0.25 + 0.3 * 1 / (1 - exp(0.04 * age)) + 0.1 * ily;
  p = p * (female ? 1.25 : 0.75);
  p = std::max(p, 0.0);
  p = std::min(p, 1.0);
  return p;
}

//' @title Use three inputs to predict response using Rcpp.
//' @description The prediction model is described in http://www.babelgraph.org/wp/?p=358.
//' @param age the first predictor (numeric)
//' @param female the second predictor (logical)
//' @param ily the third predictor (logical)
//' @return a random sample of size \code{n}
//' @examples
//' \dontrun{
//' data(data)
//' attach(data)
//' res <- vaccC(age,female,ily)
//' }
//' @export
// [[Rcpp::export]]
NumericVector vaccC(NumericVector age, LogicalVector female,
                    LogicalVector ily) {
  int n = age.size();
  NumericVector out(n);
  for(int i = 0; i < n; ++i) {
    out[i] = vacc3a(age[i], female[i], ily[i]);
  }
  return out;
}
