unit.test <- function(test,pass,fail){
	# test: logical expression to evaluate
	# pass: error to print if true
	# fail: error to print if false
	if(test){
		print(pass)
	} else{print(fail)}

}

fisher.r.to.z <- function(r){
  r <- 0.5*(log(1+r) - log(1-r))
  return(r)
}

colSD <- function(x){
  return(apply(x,2,sd))
}

fda <- function(x1,x2,nperms = 1000){
	# x1: N-by-P matrix where there are N observational units and P are features, probably time
	# x2: M-by-P matrix where there are N observational units and P are features, probably time
	# perform functional data analysis, where you compute the mean difference between column
	# means of x1 and x2 and compare them to permuted versions of x1 and x2 with the observational units
	# switched	


}

source.save <- function(script,output){
	# wrapper for source function that saves output and input of script
	file.create(output)
	con <- file(output)
	#sink(con, append=TRUE)
	sink(con, type="output")
	source(script)
	sink() 
	sink(type="output")

}

quiet <- function(x) { 
	# https://stackoverflow.com/questions/34208564/how-to-hide-or-disable-in-function-printed-message/34208658#34208658
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 