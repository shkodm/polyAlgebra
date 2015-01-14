#' Package for analitical calculations on multiple variable polynomials
#'
#' @docType package
#' @name polyAlgebra-package
#' @rdname polyAlgebra-package
NULL

#' @export
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

#' @export
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


ToC.numeric = function(x) {as.character(x)}

#' @export
PV = function(x,...,sep="") {
  if (is.character(x)) {
    x = paste(x,...,sep=sep);
  } else {
    x = x
  }
  new("gvector",vec=lapply(x,pAlg),dim=length(x))
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
	

#' @export
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

#' @export
print.pAlg = function(p)
{
	class(p) = "data.frame";
	print(p)
	print(attr(p,"var"))
}

#' @export
"+.pAlg" <- function(p1,p2){
  if (is.numeric(p1)) p1 = pAlg(p1)
  if (is.numeric(p2)) p2 = pAlg(p2)
  p = rbind(p1,p2)
	p = aggregate(p)
	p
}

#' @export
"-.pAlg" <- function(p1,p2){
  if (is.numeric(p1)) p1 = pAlg(p1)
  if (is.numeric(p2)) p2 = pAlg(p2)
  p2 = p2 * (-1)
  p = rbind(p1,p2)
  p = aggregate(p)
  p
}


#' @export
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


#' @export
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


 
#' @export
is.zero = function (x, ...) UseMethod("is.zero")

#' @export
is.zero.pAlg = function(p) all(p$.M == 0)

#' @export
is.zero.numeric = function(p) p == 0

#' @export
is.zero.gvector = function(x) gapply(x, is.zero, simplify=TRUE)

#' @export
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

#' @export
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


#' @export
subst = function (obj_, ...) UseMethod("subst")

#' @export
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


#' @export
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

