#include <Rcpp.h>

#include <sstream>

#define tostr( x ) static_cast< std::ostringstream & >( \
( std::ostringstream() << std::dec << x ) ).str()              
  
using namespace Rcpp;

// TODO: Add taking number into account (gcd for an array), when generation of sequence
// TODO: Do something with proper +/- interaction when summing
// TODO: Generally wrap it around better


// [[Rcpp::export]]
std::string dfToString(NumericMatrix coeff,List lsnames)
{
  std::string res = "";
  for (int i = 0; i < coeff.nrow(); i++)
  {
    if (coeff(i,lsnames.size() - 1) != 0) 
    {
      std::string temp = tostr(coeff(i,lsnames.size() - 1));
      
      for (int j = 0; j < (lsnames.size() - 1); j++)
      {
        std::string coff = tostr(coeff(i,j));
        switch (int(coeff(i,j))){
          case -1:
            temp = temp + "/" + as<std::string>(lsnames[j]);
            break;
          case 0: break;
          case 1:
            temp = temp + "*" + as<std::string>(lsnames[j]); 
            break;
          case 2:
            temp = temp + "*" + as<std::string>(lsnames[j]) + "*" + as<std::string>(lsnames[j]);
            break;
          default:
            temp = temp + "*" + "pow(" + as<std::string> (lsnames[j])  + "," + coff + ")";
            break;
         }
      }
      if (((coeff(i,lsnames.size() - 1)) == 1) && (temp.length() != 1) && (temp.at(1) != '/')) temp.erase(0,2);
      if (((coeff(i,lsnames.size() - 1)) == - 1) && (temp.length() != 2) && (temp.at(2) != '/')) temp.erase(1,2);
      if (res.length() == 0) res = temp;  
      else { if (temp.at(0) == '-') res = res + temp;
      else res = res + "+" + temp; }
    }
  }
  return res;
} 

  // [[Rcpp::export]]
std::string fastMult(NumericMatrix coeff,List lsnames) 
{
 // some supplementary stuff
  int mdiv_count = 0;
  int divs = -1;
  int eqzero = 0;
  int max_ind = -1;
  int temp_it = 0;
  NumericMatrix remaind = clone(coeff); 
  //we clone, since normally rcpp types are passed by pointer
  for (int i = 0; i < (lsnames.length() - 1) ; i++)
  {
      int countt = 0;
      int div_count =0;
      for (int j = 0; j < coeff.nrow(); j ++)
      {
        if ((coeff(j,i) != 0) && (coeff(j,lsnames.length() - 1) != 0 )) 
        {
          countt++;
          if ((coeff(j,i)) < 0) div_count++; 
        }
      }
      if (countt > temp_it){
        temp_it = countt;
        max_ind = i;
        mdiv_count = div_count;
      }
  }
  if ((2*mdiv_count) != (temp_it)) temp_it = - (2*mdiv_count - temp_it) / abs(2*mdiv_count - temp_it); 
  else temp_it = 1;
  int divs_ind = 0; 
  for (int j = 0; j < coeff.nrow(); j++ )
  {
    if ((coeff(j,lsnames.length() - 1) != 0 ) && (((coeff(j,max_ind) >= 1) && (temp_it == 1)) || ((coeff(j,max_ind) <= -1) && (temp_it == -1))))
    {
      coeff(j,max_ind) = coeff(j,max_ind) - temp_it;
      remaind(j,lsnames.length() - 1) = 0;
      divs_ind = divs_ind + 1; 
      eqzero++;
      divs = j;
    }
    else{
     coeff(j,lsnames.length() - 1) = 0;
     if ((remaind(j,lsnames.length() - 1)) == 0) eqzero++;
    }
  }
  if (divs_ind == 0) return dfToString(remaind,lsnames); // It should be here
  if (divs_ind == 1) {
    coeff(divs,max_ind) = coeff(divs,max_ind) + temp_it;
    std::string res = dfToString(coeff,lsnames);
    if (eqzero != coeff.nrow()) return (res + " + (" + dfToString(remaind,lsnames)+")");
    else return res;
  }
  // TODO: it is dummy, change it!
 if (temp_it < 0) {
    if (eqzero != coeff.nrow()) return  ("(" + fastMult(coeff,lsnames) + ")" + "/" + as<std::string>(lsnames[max_ind])  + "+" + fastMult(remaind,lsnames));
    else return "(" + fastMult(coeff,lsnames) + ")/" + as<std::string>(lsnames[max_ind]);
 } 
 else {
 if (eqzero != coeff.nrow()) return as<std::string>(lsnames[max_ind]) + "*" + "(" + fastMult(coeff,lsnames) + ")" + " + " + fastMult(remaind,lsnames); // <- that was working
    else return as<std::string>(lsnames[max_ind]) + "*" + "(" + fastMult(coeff,lsnames) + ")";

  }
 }
