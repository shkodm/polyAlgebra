
#' @export
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

# #' @export
# ToC.pAlg = function(p,float=TRUE, minimal=1e-10)
# {
# 	nToC(p, min=minimal,float=float)
# #	oToC(p, minimal=minimal,float=float)
# }

#' @export
ToC.pAlg = function(f,float=TRUE,minimal=1e-10)
{
  if (nrow(f) < 1) {
     			ret = "0"
     		}
  else {
    lsnames = as.list(names(f))
    matrixa = data.matrix(f)
    fastMult(matrixa,lsnames)
  }
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

M.max = 10
M.val = outer(1:M.max,1:M.max,"/")
M.PV = outer(1:M.max,1:M.max,function(x,y) paste(x,"/",y) )
#M.PV = PV(M.str)
M.sel = !duplicated(as.vector(M.val))
M.val = M.val[M.sel]
M.PV = M.PV[M.sel]


#' @export
C = function(x,y,...,eq=" = ",sep) {
  x = ToC(x,...)
  if (missing(y)) {
    if (missing(sep)) sep = ""
    cat(x,sep=sep);
  } else {
    if (missing(sep)) sep = ";\n"
    y = ToC(y,...)
    cat(paste0(x,eq,y,sep),sep="")
  }
}

#' @export
is.int = function(x,min=1e-6) {
	abs(x - round(x)) < min
}

#' @export
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

# nToC = function(tab, bracket=FALSE,min=1e-6, second=FALSE, float=TRUE) {
# 	tab = tab[abs(tab$.M) > min,,drop=FALSE]
# 	if (nrow(tab) < 1) {
# 		if (second) {
# 			ret = " + 0"
# 		} else {
# 			ret = "0"
# 		}
# 	} else {
# 		tab = tab[order(tab$.M,decreasing=TRUE),,drop=FALSE]
# 		i1=colSums(tab > 0)
# 		i2=colSums(tab < 0)
# 		Md = data.frame(
# 			val = c(1:36,1/(1:36)),
# 			str = {
# 				if (float) {
# 					str = paste(c(1:36,1:36),rep(c(".","."),each=36),sep="")
# 				} else {
# 					str = paste(c(1:36,1:36),rep(c("","."),each=36),sep="")
# 				}
# 			},
# 			positive = rep(c(TRUE,FALSE),each=36)
# 		)
# 		Md = Md[c(36:1,1:36+36),]
# 		Md.val = tab$.M
# 		Md = rbind(Md, data.frame(
# 			val = Md.val,
# 			str = paste(tab$.M,"",sep=""),
# 			positive = TRUE
# 		))
# 		Md = no.ones(Md)
# 		if (nrow(Md) > 0) {
# 			i3t = divisible(tab$.M, Md$val)
# 			i3 = colSums(i3t)
# 		} else {
# 			i3 = 0
# 		}
# 		i1[".M"] = -1
# 		i2[".M"] = -1
# 		if (any(c(i1,i2,i3) > 0)) {
# 			wh = which.max(c(max(i3),max(i2),max(i1)))
# 			
# 			if (wh == 1) {
# 				i = which.max(i3)
# 				sel = i3t[,i]
# 				positive = Md$positive[i]
# 				ntab = tab[sel,,drop=FALSE]
# 				ntab$.M = ntab$.M / Md$val[i]
# 				pull = Md$str[i]
# 			} else if (wh == 3) {
# 				i = which.max(i1)
# 				sel = tab[,i] > 0
# 				positive=TRUE
# 				ntab = tab[sel,,drop=FALSE]
# 				ntab[,i] = ntab[,i]-1
# 				pull = names(tab)[i]
# 			} else if (wh == 2) {
# 				i = which.max(i2)
# 				sel = tab[,i] < 0
# 				positive=FALSE
# 				ntab = tab[sel,,drop=FALSE]
# 				ntab[,i] = ntab[,i]+1
# 				pull = names(tab)[i]
# 			}
# 			if (any(!sel)) {
# 				v1 = nToC(ntab,bracket=T,second=TRUE,float=float)
# 			} else {
# 				v1 = nToC(ntab,bracket=T,second=second,float=float)
# 			}
# 			if (positive) {
# 				if (v1 == "1") {
# 					v1 = paste(pull,sep="")
# 				} else if (v1 == " + 1") {
# 					v1 = paste(" + ",pull,sep="")
# 				} else if (v1 == "-1") {
# 					v1 = paste("-",pull,sep="")
# 				} else if (v1 == " - 1") {
# 					v1 = paste(" - ",pull,sep="")
# 				} else {
# 					v1 = paste(v1,"*",pull,sep="")
# 				}
# 			} else {
# 				v1 = paste(v1,"/",pull,sep="")
# 			}
# 			if (any(!sel)) {
# 				if (bracket) {
# 					v2 = nToC(tab[!sel,,drop=FALSE],second=FALSE,float=float)
# 					if (second) {
# 						ret = paste(" + ( ",v2,v1," )",sep="")
# 					} else {
# 						ret = paste("( ",v2,v1," )",sep="")
# 					}
# 				} else {
# 					v2 = nToC(tab[!sel,,drop=FALSE],second=second,float=float)
# 					ret = paste(v2,v1,sep="")
# 				}
# 			} else {
# 				ret = v1
# 			}
# 		} else {
# 			v = sum(tab$.M)
# 			if (abs(round(v) - v) < min) {
# 				v = round(v)
# 				ret = sprintf("%d",abs(v))
# 			} else {
# 				ret = sprintf("%.16f",abs(v))
# 			}
# 			if (second) {
# 				if (v < 0) {
# 					ret = paste(" - ",ret,sep="")
# 				} else {
# 					ret = paste(" + ",ret,sep="")
# 				}
# 			} else {
# 				if (v < 0) {
# 					ret = paste("-",ret,sep="")
# 				} 
# 			}
# 		}
# 	}
# 	ret
# }


 
#' @export
ToC.gvector = function(x,...) 
{ gapply(x, ToC, ..., simplify=TRUE) }
#{
 # lsnames = as.list(names(x[[1]]))
#  matrixa = data.matrix(x[[1]])
#  fastMult(matrixa,lsnames)
#}

#' @export
ToC.numeric = function(x) {as.character(x)}



