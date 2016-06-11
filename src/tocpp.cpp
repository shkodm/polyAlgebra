#include <Rcpp.h>

#include <sstream>

#include <math.h>

#define tostr( x ) ((is_int(x))? (static_cast< std::ostringstream & >(\
( std::ostringstream() << std::dec << x ) ).str()): (static_cast< std::ostringstream & >( \
    ( std::ostringstream() << std::dec << std::setprecision(16) << x ) ).str()))

#define is_int( x ) (fabs(x - round(x)) < 1e-10)

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix ugly_order(NumericMatrix x, Function f)
{
  NumericMatrix res = f(x);
  return res;
}
//TODO: add check is decimal is finite, so unnecessary simplification to fraction can be avoided

// [[Rcpp::export]]
std::string dfToString(NumericMatrix coeff_in,List lsnames,bool flt, Function foo)
{
  std::string res = "";
  std::string ddot = "";
  std::string ssign;
  NumericMatrix coeff = ugly_order(coeff_in,foo);
  int ncol = coeff.ncol();
  int nrow = coeff.nrow();
  if (as<std::string>(lsnames[0]) == ".M") return tostr(coeff(0,0));
  for (int i = 0; i < nrow; i ++)
  {
    if (coeff(i,ncol - 1) != 0) 
    {
      std::string temp = "";
      std::string ssign = "*";
      if ((is_int(coeff(i,ncol - 1))) && (flt)) ddot = ".";
      else ddot = "";
      for (int j = 0; j < (ncol - 1); j++)
      {
        std::string coff = tostr(coeff(i,j));
        if ((!temp.empty()) && (temp.at(0) == '/')) ssign = ""; //it is here, when it is empty
        else ssign = "*";
        if (is_int(coeff(i,j))){   
          switch (int(coeff(i,j))){
          case -1:
            temp ="/" +  as<std::string>(lsnames[j]) + ssign + temp;
            break;
          case 0:
            break;
          case 1:
            temp = as<std::string>(lsnames[j]) + ssign + temp;
            break;
          case 2:
            temp = as<std::string>(lsnames[j]) + "*" + as<std::string>(lsnames[j]) + ssign + temp;
            break;
          default:
            temp = "pow(" + as<std::string> (lsnames[j])  + "," + coff + ")" + ssign + temp;
          break;
          }}
        else temp = "pow(" + as<std::string> (lsnames[j])  + "," + coff + ")" + ssign + temp;
      }
      if (tostr(fabs(coeff(i,ncol - 1))) == "1") {
        if ((!temp.empty()) && (temp.at(temp.size() - 1) == '*')) temp.erase(temp.end() - 1);
        if (temp.empty()) temp = tostr(coeff(i,ncol - 1)); 
        else if (temp.at(0) == '/') temp = tostr(coeff(i, ncol - 1)) + ddot + temp;
        else if ((coeff(i,ncol - 1 ) < 0)) temp = "-" + temp;
      }
      else {
        if ((!temp.empty()) && (temp.at(0) == '/')) 
        {
          temp.erase(temp.end() - 1);
          temp = tostr(coeff(i, ncol - 1)) + ddot + temp;
        }
        else if (coeff(i,ncol - 1) < 0) temp ="-" + temp + tostr(fabs(coeff(i,ncol - 1))) + ddot;
        else temp = temp + tostr(coeff(i, ncol - 1)) + ddot;
      }
      if (((temp == "1") || (temp == "1.") || (temp == "-1.") || (temp == "-1")) && (res.length() != 0)) 
      {
        if (res.at(0) == '-') 
        {
          res.insert(0," ");
          res.insert(2," ");
          res.insert(0,temp);
        }
        else res.insert(0,(temp + " + "));
      } else { if (res.empty()) res = temp;  
      else { 
        if (temp.at(0) == '-') {
          temp.insert(1," ");
          res = res + " " + temp;
        }
        else res = res + " + " + temp; }
      } 
    } }
  return res;
  
}

// [[Rcpp::export]]
std::string fastMult(NumericMatrix coeff,List lsnames,bool flt,Function foo) 
{
  NumericMatrix remaind = clone(coeff);
  
  int mdiv_count = 0;
  double divs = -1; 
  double ddivs = -1;
  int eqzero = 0;
  int max_ind = 0;
  int temp_it = 0;
  int qmax_ind = 0;
  int dmax_ind = 0;
  double curr_el; 
  std::string sgn;
  int ncol = coeff.ncol();
  int nrow = coeff.nrow();
  std::string l_bracket = "";
  std::string r_bracket = "";
  if ( (as<std::string>(lsnames(ncol - 1)) == ".o")) 
  {
    l_bracket = l_bracket + "( ";
    r_bracket = r_bracket + " )";
    lsnames(ncol - 1) = ".M";
  }
  
  for (int divsor = 2; divsor < 37; divsor++ )
  {
    temp_it = 0;
    mdiv_count = 0;
    for (int j = 0 ; j < nrow; j ++) 
    {
      double curr_el = coeff(j, ncol - 1);
      int check_it = int(coeff(j, ncol - 1));
      if (is_int(curr_el)) 
      {
        
        if ((int(curr_el) != 0) && (is_int(curr_el/divsor))) temp_it++; 
      }
      else 
      {
        if (is_int(curr_el*divsor)) mdiv_count++;
      }
    }
    if ((temp_it >= dmax_ind) && (temp_it > 0)) 
    {
      dmax_ind = temp_it;
      divs = divsor;
    }
    if ((mdiv_count > qmax_ind) && (mdiv_count > 0))
    {
      qmax_ind = mdiv_count;
      ddivs = divsor;
    }
  }
  
  double rdivs = -1;
  int rmax_ind = -1;
  for (int i = 0; i < nrow; i++)
  {
       temp_it = 0;
       if ((std::fabs(coeff(i, ncol - 1 )) != 1.) && (coeff(i, ncol - 1 ) != 0)) {
         for (int j = 0; j < nrow; j++) {
           if ((coeff(j, ncol - 1 ) != 0) && (is_int(coeff(j, ncol - 1)/coeff(i, ncol - 1))) && (is_int(coeff(j, ncol - 1)) ? (is_int(coeff(i,ncol - 1))): (!is_int(coeff(i, ncol - 1))))) temp_it++;} }
       if ((temp_it >= rmax_ind) && (temp_it > 0))
       {
         rdivs = coeff(i, ncol - 1);
         rmax_ind = temp_it;
       }
  }
     
  
  mdiv_count = 0;
  temp_it = 0;
  for (int i = 0; i < (ncol - 1) ; i++)
  {
    int countt = 0;
    int div_count = 0;
    for (int j = 0; j < nrow; j ++)
    {
      if (double(coeff(j,i) != 0.) && (double(coeff(j,ncol - 1)) != 0. )) 
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
  
 bool self_div = ((rmax_ind > 1) && (rmax_ind > qmax_ind) && (rmax_ind >= dmax_ind));
///self_div = false; // if something is slightly off
  if (qmax_ind >= dmax_ind) divs = ddivs;
  if (self_div)
  {
    dmax_ind = rmax_ind;
    divs = std::fabs(rdivs);
  }
  if ((qmax_ind >= temp_it) || (dmax_ind >= temp_it)) {
    if ((std::fabs(divs) != 1) && ((dmax_ind != 1) || (qmax_ind > 0))) {
      for (int j = 0; j < nrow; j++)
      {
        double curr_el = coeff(j, ncol - 1);
        if ((curr_el != 0.) && ((((is_int(curr_el/divs))) && (dmax_ind > qmax_ind) && (is_int(curr_el)) && (!self_div))  // <- i'm not proud of this
                                  ||  (((is_int(curr_el*divs) ) && (qmax_ind >= dmax_ind) 
                                          && (!(is_int(curr_el))) && (!self_div))) || ((self_div) && (is_int(curr_el/divs))))) 
        {
          if ((dmax_ind > qmax_ind ) || (self_div)) coeff(j, ncol - 1) = coeff(j, ncol - 1) / divs;
          else coeff(j, ncol - 1) = coeff(j, ncol - 1) * divs;
          remaind(j,ncol - 1) = 0;
          eqzero++;
        }
        else  {     
          coeff(j, ncol - 1) = 0;
          if ((remaind(j,ncol - 1)) == 0) eqzero++;
        }
      }
      std::string d_dot ="";
      if ((is_int(divs)) && (flt)) d_dot = ".";
      if ((dmax_ind > qmax_ind) || (self_div)) sgn = "*";
      else sgn = "/";
      if (qmax_ind != 1) lsnames(ncol - 1) = ".o";
      if (eqzero == nrow - 1) 
      {
        std::string  singleDivs = fastMult(coeff,lsnames,flt,foo);
        if (singleDivs.at(0) == '-') {
          singleDivs.insert(1," ");
          return (l_bracket + fastMult(remaind,lsnames,flt,foo) + " " + singleDivs + sgn +  tostr(divs) + d_dot + r_bracket);
        }
        else return (l_bracket + fastMult(remaind,lsnames,flt,foo) + " + " + singleDivs + sgn +  tostr(divs) + d_dot + r_bracket); 
      }
      else if (eqzero != nrow)  return (l_bracket + fastMult(remaind,lsnames,flt,foo) + " + " + fastMult(coeff,lsnames,flt,foo) + sgn +  tostr(divs) + d_dot + r_bracket); 
      else return (fastMult(coeff,lsnames,flt,foo) + sgn + tostr(divs) + d_dot);
    }
    else return (l_bracket + dfToString(coeff,lsnames,flt,foo) + r_bracket); 
  }
  else {
    
    eqzero = 0;
    divs = -1;
    
    if (temp_it == 1)  return (l_bracket + dfToString(coeff,lsnames,flt,foo) + r_bracket); 
    if ((2*mdiv_count) != (temp_it)) temp_it = - (2*mdiv_count - temp_it) / abs(2*mdiv_count - temp_it); 
    else temp_it = 1;
    int divs_ind = 0; 
    for (int j = 0; j < nrow; j++ )
    {
      if (double(coeff(j,ncol - 1) != 0. ) && (((coeff(j,max_ind) >= 1) && (temp_it == 1)) || ((coeff(j,max_ind) <= -1) && (temp_it == -1))))
      {
        coeff(j,max_ind) = coeff(j,max_ind) - temp_it;
        remaind(j,ncol - 1) = 0;
        divs_ind = divs_ind + 1; 
        eqzero++;
        divs = j;
      }
      else{
        coeff(j,ncol - 1) = 0;
        if ((remaind(j,ncol - 1)) == 0) eqzero++;
      }
    }
    if (as<std::string>(lsnames[max_ind]) == ".M") return (l_bracket + dfToString(coeff,lsnames,flt,foo) + r_bracket);
    if (divs_ind == 0) return  (l_bracket + dfToString(remaind,lsnames,flt,foo) + r_bracket);  
    if (temp_it < 0) sgn = "/";
    else sgn = "*";
    lsnames(ncol - 1) = ".o";
    if (eqzero == nrow - 1) 
    {
      std::string  singleDivs = fastMult(coeff,lsnames,flt,foo);
      if (singleDivs.at(0) == '-') {
        singleDivs.insert(1," ");
        return (l_bracket + fastMult(remaind,lsnames,flt,foo) + " " + singleDivs + sgn +  as<std::string>(lsnames[max_ind]) + r_bracket);
      }
      else return (l_bracket + fastMult(remaind,lsnames,flt,foo) + " + " + singleDivs + sgn +  as<std::string>(lsnames[max_ind]) + r_bracket); 
    }
    else if (eqzero != nrow) return  ( l_bracket + fastMult(remaind,lsnames,flt,foo) +  " + " + fastMult(coeff,lsnames,flt,foo)  + sgn + as<std::string>(lsnames[max_ind]) + r_bracket);
    else return  (fastMult(coeff,lsnames,flt,foo) + sgn + as<std::string>(lsnames[max_ind])); 
  } 
}
