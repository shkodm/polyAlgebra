#include <Rcpp.h>
#include <algorithm>
#include <sstream>
#include <string>
using namespace Rcpp;
#define is_int( x ) (fabs(x - round(x)) < 1e-10)
#define tostr( x ) ((is_int(x))? (static_cast< std::ostringstream & >(                    \
( std::ostringstream() << std::dec << x ) ).str()): (static_cast< std::ostringstream & >( \
    ( std::ostringstream() << std::dec << std::setprecision(16) << x ) ).str()))

// [[Rcpp::export]]
NumericVector poly_aggregate(NumericMatrix a, NumericVector dotM)
{
  std::vector <std::string> val_buffer;
  int ncol = a.ncol();
  int nrow = a.nrow();
  for (int i = 0; i < nrow; ++i)
  {
    std::string curr;
    for (int j = 0; j < ncol; ++j) curr = curr + tostr(a(i,j)); 
    std::vector <std::string>::iterator pos;
    if (!val_buffer.empty()) pos = (std::find(val_buffer.begin(), val_buffer.end(), curr));
    if (val_buffer.empty() || (pos == val_buffer.end())) val_buffer.push_back(curr);
    else {
      dotM[pos - val_buffer.begin()] = dotM[pos - val_buffer.begin()] + dotM[i];  
      dotM[i] = 0;
    }
  }
 return dotM;
}