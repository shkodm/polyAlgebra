#include <Rcpp.h>

#include <sstream>

#define tostr( x ) static_cast< std::ostringstream & >( \
( std::ostringstream() << std::dec << x ) ).str()              

using namespace Rcpp;


// [[Rcpp::export]]
std::string dfToString(NumericMatrix coeff,List lsnames)
{
  std::string res = "";
  if (as<std::string>(lsnames[0]) == ".M") return tostr(coeff(0,0));
  for (int i = 0; i < coeff.nrow(); i ++)
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
        case 0:
          break;
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
      else { if (temp.at(0) == '-') {
        temp.insert(1," ");
        res = res + " " + temp;
      }
      else res = res + " + " + temp; }
    }
  }
  return res;
} 


// [[Rcpp::export]]
std::string fastMult(NumericMatrix coeff,List lsnames) 
{
  NumericMatrix remaind = clone(coeff);
  
  int mdiv_count = 0;
  int divs = -1;
  int ddivs = -1;
  int eqzero = 0;
  int max_ind = 0;
  int temp_it = 0;
  int qmax_ind = 0;
  int dmax_ind = 0;
  std::string sgn;
  
  std::string l_bracket = "";
  std::string r_bracket = "";
  if ( (as<std::string>(lsnames(lsnames.length() - 1)) == ".o")) 
  {
    l_bracket = l_bracket + "( ";
    r_bracket = r_bracket + " )";
    lsnames(lsnames.length() - 1) = ".M";
   }
  
  for (int divsor = 2; divsor < 37; divsor ++)
  {
    temp_it = 0;
    mdiv_count = 0;
    for (int j = 0 ; j < coeff.nrow(); j ++)
      if (double(int(coeff(j, lsnames.length() - 1))) == coeff(j, lsnames.length() - 1)) {
        if ((int(coeff(j, lsnames.length() - 1)) != 0) && (int((coeff(j, lsnames.length() - 1))) % divsor == 0)) temp_it++; }
        else 
        {
          if ((((coeff(j, lsnames.length() - 1))*divsor) == (double(int((coeff(j, lsnames.length() - 1))*divsor))))) mdiv_count++;
        }
        if (temp_it > dmax_ind) 
        {
          dmax_ind = temp_it;
          divs = divsor;
        }
        if (mdiv_count > qmax_ind)
        {
         qmax_ind = mdiv_count;
         ddivs = divsor;
        }
  }
  
  mdiv_count = 0;
  temp_it = 0;
  for (int i = 0; i < (lsnames.length() - 1) ; i++)
  {
    int countt = 0;
    int div_count = 0;
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
  
  if (qmax_ind >= dmax_ind) divs = ddivs;
  if ((qmax_ind >= temp_it) || (dmax_ind >= temp_it)) {
    if ((divs != -1) && ((dmax_ind != 1) || (qmax_ind > 0))) {
      for (int j = 0; j < coeff.nrow(); j++)
      {
        if ((double(coeff(j, lsnames.length() - 1)) != 0.) && ((((int((coeff(j, lsnames.length() - 1))) % divs == 0)) && (dmax_ind > qmax_ind))  
                ||  ((((coeff(j, lsnames.length() - 1))*divs) == double(int((coeff(j, lsnames.length() - 1))*divs))) && (qmax_ind >= dmax_ind) 
                       && (double(int(coeff(j, lsnames.length() - 1))) != coeff(j, lsnames.length() - 1)) ))) 
        {
          if (dmax_ind > qmax_ind ) coeff(j, lsnames.length() - 1) = coeff(j, lsnames.length() - 1) / divs;
          else coeff(j, lsnames.length() - 1) = coeff(j, lsnames.length() - 1) * divs;
          remaind(j,lsnames.length() - 1) = 0;
          eqzero++;
        }
        else  {     
          coeff(j, lsnames.length() - 1) = 0;
          if ((remaind(j,lsnames.length() - 1)) == 0) eqzero++;
        }
      }
    std::string d_dot =".";
    if (divs == (double(int(divs)))) d_dot = ".";
    if (dmax_ind > qmax_ind) sgn = "*";
    else sgn = "/";
    if (qmax_ind != 1) lsnames(lsnames.length() - 1) = ".o";
    if (eqzero != coeff.nrow())  return (l_bracket + fastMult(remaind,lsnames) + " + " + fastMult(coeff,lsnames) + sgn +  tostr(divs) + d_dot + r_bracket); 
    else return (fastMult(coeff,lsnames) + sgn + tostr(divs) + d_dot);
    }
    else return (l_bracket + dfToString(coeff,lsnames) + r_bracket);
  }
  else {
    
    eqzero = 0;
    divs = -1;
    
    if (temp_it == 1) return (l_bracket + dfToString(coeff,lsnames) + r_bracket);
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
    if (as<std::string>(lsnames[max_ind]) == ".M") return (l_bracket + dfToString(coeff,lsnames) + r_bracket);
    if (divs_ind == 0) return  (l_bracket + dfToString(remaind,lsnames) + r_bracket); 
    if (temp_it < 0) sgn = "/";
    else sgn = "*";
    lsnames(lsnames.length() - 1) = ".o";
    if (eqzero != coeff.nrow()) return  ( l_bracket + fastMult(remaind,lsnames) +  " + " + fastMult(coeff,lsnames)  + sgn + as<std::string>(lsnames[max_ind]) + r_bracket);
    else return  (fastMult(coeff,lsnames) + sgn + as<std::string>(lsnames[max_ind])); 
  } 
}
