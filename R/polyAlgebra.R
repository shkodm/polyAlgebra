#' Package for analitical calculations on multiple variable polynomials
#'
#' @docType package
#' @name polyAlgebra-package
#' @rdname polyAlgebra-package
#' @useDynLib polyAlgebra
NULL

#' Class that represents a polynomial function
#' 
#' @slot tab Array keepAlging all the information about the function (see below)
#' @exportClass pAlg
pAlg = setClass("pAlg", representation(tab="array"))


pAlg = function(a){
	if (is.character(a)) {
		ret = data.frame(.M = 1)
		ret[,a] = 1;
		attr(ret,"var") = a
	} else if (is.numeric(a)) {
		ret = data.frame(.M = a)
#		if (a == 0) { ret = data.frame(.M = c()) }
		attr(ret,"var") = c()
	} else {
		stop("Unknow type in pAlg()\n");
	}
	finish.pAlg(ret)
}	

scalar = function(v,a){
	
	if ("pAlg" %in% class(a)) { a = as.factor(attr(a,"var")) }
	if (is.factor(a)) {
		ret = matrix(0, 1, nlevels(a)+1)
		ret = data.frame(ret)
		names(ret) = c(".M", levels(a))
		ret$.M = v;
		attr(ret,"var") = levels(a)
	} else {stop("a have to be a factor in scalar()");}
	finish.pAlg(ret)
}	


aggregate.pAlg = function(p)
{
	if (nrow(p) > 0) {
                class(p) = "data.frame"
		if (length(attr(p,"var")) == 0) {
	                ret = p[1,,drop=FALSE]
			ret$.M = sum(p$.M)
		} else {
	                ret = aggregate(p[,".M",drop=FALSE], p[,attr(p,"var"),drop=FALSE], sum)
		}
                ret = ret[ret$.M != 0,,drop=FALSE]
                attr(ret,"var") = attr(p,"var")
                p = finish.pAlg(ret)
	}
	if (nrow(p) == 0) {
		p = pAlg(0)
	}
	p
}

finish.pAlg = function(p) {
	v = names(p)
	sel = v %in% ".M"
	if (all(!sel)) stop("There should be a .M in pAlg object")
	class(p) = c("pAlg","data.frame")
	attr(p,"var") = v[!sel]
	p
}
	

rbind.pAlg = function(p1,p2)
{
	if (is.null(p1)) return(p2)
	if (is.null(p2)) return(p1)
	if (nrow(p1) == 0) return(p2);
	if (nrow(p2) == 0) return(p1);
	class(p1) = "data.frame"
	class(p2) = "data.frame"
	col = names(p1)[!names(p1) %in% names(p2)]
	p2[,col] = 0;
	col = names(p2)[!names(p2) %in% names(p1)]
	p1[,col] = 0;
	ret = rbind(p1,p2)
	attr(ret,"var") = c(attr(p1,"var"),col)
	finish.pAlg(ret)
}

print.pAlg = function(p)
{
	class(p) = "data.frame";
	print(p)
	print(attr(p,"var"))
}

"+.pAlg" <- function(p1,p2){
  if (is.numeric(p1)) p1 = pAlg(p1)
  if (is.numeric(p2)) p2 = pAlg(p2)
  p = rbind(p1,p2)
	p = aggregate(p)
	p
}

"^.pAlg" <- function(p1,p2){
	if (is.numeric(p2)) {
		p = p1;
		for (i in names(p1)) {
			if (i == ".M") {
				p[,i] = p[,i] ^ p2
			} else {
				p[,i] = p[,i] * p2
			}
		}
	} else {
		stop("non numeric power in ^.pAlg")
	}
	p = aggregate(p)
	p
}


"*.pAlg" <- function(p1,p2){
#	cat("-- *.pAlg -----\n");
  if (is.numeric(p1)) {
    tmp = p1
    p1 = p2
    p2 = tmp
  }
	if (! "pAlg" %in%  class(p2)) {
		if (is.numeric(p2)) {
			if (length(p2) > 1) stop("pAlg only multiply by numeric length 1\n")
			p1$.M = p1$.M * p2
			p1 = aggregate(p1)
			return(p1)
		}
		stop("Unknown type in pAlg multyply\n");
	}
#	cat("----------------------------------------------------------------\n");
#	print(c(nrow(p1),nrow(p2)))
#	print(p1)
#	print(p2)

	if ((nrow(p1) != 0)&&(nrow(p2) != 0)) {               
                i = rep(1:nrow(p1),each =nrow(p2))
                j = rep(1:nrow(p2),times=nrow(p1))
                class(p1) = "data.frame"
                class(p2) = "data.frame"
		col = names(p1)[!names(p1) %in% names(p2)]
		p2[,col] = 0;
		col = names(p2)[!names(p2) %in% names(p1)]
		p1[,col] = 0;
                v = union(names(p1), names(p2))
		v = setdiff(v , ".M")
		if (length(v) != 0) {
	                p = p1[i,v,drop=FALSE] + p2[j,v,drop=FALSE]
        	        p$.M = p1$.M[i] * p2$.M[j]
		} else {
			p = pAlg(p1$.M[i] * p2$.M[j])
		}
                attr(p,"var") = v
                class(p) = c("pAlg","data.frame")
                p = aggregate(p)
	} else {
		if (nrow(p1) == 0) p = p1;
		if (nrow(p2) == 0) p = p2;
	}	
	p
}


 
ToC = function (x, ...) 
UseMethod("ToC")


ToC_row = function(x,float=TRUE,minimal=1e-10)
{
	if (abs(x[".M"]) < minimal) x[".M"] = 0;
	if (x[".M"] != 0) {
		val = x[".M"]
		x[".M"] = abs(x[".M"])
		if (x[".M"] == 1) {
			ret = NULL
		} else {
			if (float) {
				ret = sprintf("%.10e",x[".M"])
			} else {
				ret = sprintf("%d",x[".M"])
			}	
#			ret = as.character(x[".M"])
		}
		v = names(x)
		for (i in 1:length(x))
		{
			vv = v[i]
			if (vv != ".M") {
				h = "";
				if (x[i] != 0) h = paste("pow(",vv,",",x[i],")")					
				if (x[i] == 1) h = vv					
				if (x[i] == 2) h = paste("(",vv,"*",vv,")",sep="")					
				if (x[i] == -1) h = paste("(1/",vv,")",sep="")					
				if (h != "") {
					if (is.null(ret)) {ret = h;} else {
						 ret = paste(ret,h,sep="*");
					}
				}
			}
		}
		if (is.null(ret)) {ret = "1"}
		ret = paste(ifelse(val>0," + "," - "), ret, sep="")
	} else {
		ret =""
	}
	ret
}

ToC.pAlg = function(p,float=TRUE, minimal=1e-10)
{
	nToC(p, min=minimal,float=float)
#	oToC(p, minimal=minimal,float=float)
}

oToC = function(p,float=TRUE, minimal=1e-10)
{
	if (nrow(p) > 0) {
		ret = apply(p,1,function(x) {ToC_row(x,float=float,minimal=minimal)})
	} else { ret = "   0"; }
	ret = paste(ret,collapse="");
	if (substr(ret,2,2) == "+") substr(ret,2,2) = " ";
	if (ret == "") { ret = "   0"; }
	ret
}

is.zero.pAlg = function(p) all(p$.M == 0)
is.zero = function (x, ...) UseMethod("is.zero")

der = function (x, ...) UseMethod("der")

der_row = function(x)
{
	val = x[".M"]
	v = names(x)
#	print(v)
	ret = NULL
	for (i in 1:length(x))
	{
		vv = v[i]
		if (vv != ".M") {
		if (x[i] > 0) {
			np = x;
			np[".M"] = np[".M"] * np[i]
			np[i] = np[i] - 1;
			np[der(vv)] = 1
#			print(np)
			ret = rbind(ret,np)
		}
		}
	}
#	print(ret)
	data.frame(ret)
}

der.pAlg = function(p)
{
	class(p) = "data.frame"
	v = attr(p,"var")
	tmp = matrix(0, nrow(p), length(v))
	tmp = data.frame(tmp)
	names(tmp) = der(v)
	v = c(v,names(tmp))
	p = cbind(p,tmp)

	ret = apply(p,1,der_row)
#	print(ret)
	ret = do.call(rbind, ret)
#	print(ret)
	attr(ret,"var")=v
	class(ret) = c("pAlg","data.frame")
#	print(ret)
        ret = aggregate(ret)
#	print(ret)
	ret
}

der.character = function (x) {
	nx = sub("\\[","_d[",x)
	nx = ifelse( x == nx, paste(x,"_d",sep=""), nx)
	nx
}



M.max = 10
M.val = outer(1:M.max,1:M.max,"/")
M.PV = outer(1:M.max,1:M.max,function(x,y) paste(x,"/",y) )
#M.PV = PV(M.str)
M.sel = !duplicated(as.vector(M.val))
M.val = M.val[M.sel]
M.PV = M.PV[M.sel]


C = function(x,...) {cat(ToC(x,...), sep="");}

is.int = function(x,min=1e-6) {
	abs(x - round(x)) < min
}

divisible = function(x,y,min=1e-6) {
	M.w = outer(x, y, "/")
	M.h = outer(!is.int(x),is.int(y),"|")
	is.int(M.w) & M.h
}

no.ones = function(tab,min=1e-6) {
	x = tab$val
	sel = pmin(abs(x - 1),abs(x),abs(x+1)) < min
	tab[!sel,,drop=FALSE]
}

nToC = function(tab, bracket=FALSE,min=1e-6, second=FALSE, float=TRUE) {
	tab = tab[abs(tab$.M) > min,,drop=FALSE]
	if (nrow(tab) < 1) {
		if (second) {
			ret = " + 0"
		} else {
			ret = "0"
		}
	} else {
		tab = tab[order(tab$.M,decreasing=TRUE),,drop=FALSE]
		i1=colSums(tab > 0)
		i2=colSums(tab < 0)
		Md = data.frame(
			val = c(1:36,1/(1:36)),
			str = {
				if (float) {
					str = paste(c(1:36,1:36),rep(c(".","."),each=36),sep="")
				} else {
					str = paste(c(1:36,1:36),rep(c("","."),each=36),sep="")
				}
			},
			positive = rep(c(TRUE,FALSE),each=36)
		)
		Md = Md[c(36:1,1:36+36),]
		Md.val = tab$.M
		Md = rbind(Md, data.frame(
			val = Md.val,
			str = paste(tab$.M,"",sep=""),
			positive = TRUE
		))
		Md = no.ones(Md)
		if (nrow(Md) > 0) {
			i3t = divisible(tab$.M, Md$val)
			i3 = colSums(i3t)
		} else {
			i3 = 0
		}
		i1[".M"] = -1
		i2[".M"] = -1
		if (any(c(i1,i2,i3) > 0)) {
			wh = which.max(c(max(i3),max(i2),max(i1)))
			
			if (wh == 1) {
				i = which.max(i3)
				sel = i3t[,i]
				positive = Md$positive[i]
				ntab = tab[sel,,drop=FALSE]
				ntab$.M = ntab$.M / Md$val[i]
				pull = Md$str[i]
			} else if (wh == 3) {
				i = which.max(i1)
				sel = tab[,i] > 0
				positive=TRUE
				ntab = tab[sel,,drop=FALSE]
				ntab[,i] = ntab[,i]-1
				pull = names(tab)[i]
			} else if (wh == 2) {
				i = which.max(i2)
				sel = tab[,i] < 0
				positive=FALSE
				ntab = tab[sel,,drop=FALSE]
				ntab[,i] = ntab[,i]+1
				pull = names(tab)[i]
			}
			if (any(!sel)) {
				v1 = nToC(ntab,bracket=T,second=TRUE,float=float)
			} else {
				v1 = nToC(ntab,bracket=T,second=second,float=float)
			}
			if (positive) {
				if (v1 == "1") {
					v1 = paste(pull,sep="")
				} else if (v1 == " + 1") {
					v1 = paste(" + ",pull,sep="")
				} else if (v1 == "-1") {
					v1 = paste("-",pull,sep="")
				} else if (v1 == " - 1") {
					v1 = paste(" - ",pull,sep="")
				} else {
					v1 = paste(v1,"*",pull,sep="")
				}
			} else {
				v1 = paste(v1,"/",pull,sep="")
			}
			if (any(!sel)) {
				if (bracket) {
					v2 = nToC(tab[!sel,,drop=FALSE],second=FALSE,float=float)
					if (second) {
						ret = paste(" + ( ",v2,v1," )",sep="")
					} else {
						ret = paste("( ",v2,v1," )",sep="")
					}
				} else {
					v2 = nToC(tab[!sel,,drop=FALSE],second=second,float=float)
					ret = paste(v2,v1,sep="")
				}
			} else {
				ret = v1
			}
		} else {
			v = sum(tab$.M)
			if (abs(round(v) - v) < min) {
				v = round(v)
				ret = sprintf("%d",abs(v))
			} else {
				ret = sprintf("%.16f",abs(v))
			}
			if (second) {
				if (v < 0) {
					ret = paste(" - ",ret,sep="")
				} else {
					ret = paste(" + ",ret,sep="")
				}
			} else {
				if (v < 0) {
					ret = paste("-",ret,sep="")
				} 
			}
		}
	}
	ret
}

subst = function (obj_, ...) UseMethod("subst")

subst.pAlg = function(obj_, ...) {
	arg = list(...)
	if (length(arg) == 0) return(obj_)
	if (is.null(names(arg))) names(arg) = rep("", length(arg))
	sel = names(arg) == ""
	narg = arg[!sel]
	for (l in arg[sel]) {
		narg = c(narg,l)
	}
	arg=narg
	if (any(names(arg) == "")) stop("All arguments to subst have to be named")
	sel = names(arg) %in% names(obj_)
	arg = arg[sel]
	if (length(arg) == 0) return(obj_)
	for (n in names(arg)) {
		v = arg[[n]]
		if (is.numeric(v)) v = pAlg(v)
		if (!"pAlg" %in% class(v)) stop("Substitutions have to be numeric of pAlg type in subst")
		arg[[n]] = v
	}
	sel = names(obj_) %in% names(arg)
	sum = pAlg(0)
	for (i in 1:nrow(obj_)) {
		K = as.matrix(obj_[i,names(arg)])
		ret = finish.pAlg(obj_[i,!sel,drop=FALSE])
		for (j in 1:length(arg)) {
			if (K[j] < 0) stop("Negative powers not supported in subst")
			if (K[j] > 0) for (l in 1:K[j]) ret = ret * arg[[j]]
		}
		sum = sum + ret
	}
	sum
}


deriv.pAlg = function(obj_, what) {
	if (!is.character(what)) stop("Can only deriv pAlg and PV with respect to variables given as strings")
	if (what %in% names(obj_)) {
		obj_$.M = obj_$.M * obj_[,what];
		obj_[,what] = obj_[,what] - 1;
	} else {
		obj_$.M = 0;
	}
	aggregate(obj_)	
}


