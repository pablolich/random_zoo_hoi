roots2 = function(a, b, c){
  x1=-b+sqrt(b^2-4*a*c)
  x2 = -b-sqrt(b^2-4*a*c)
  return(cbind(x1, x2))
}

countevents = function(data, event1, event2){
  nboth = 0
  for (i in 1:nrow(data)){
    if (data[i,1]==event1){
      if (data[i,2]==event2){
        nboth = nboth + 1
      }
    else (nboth = nboth)
    }

  }
  return(nboth/nrow(data))
}

n = 1000
a = rnorm(n, 0, 1)
b = rnorm(n, 0, 1)
c = rnorm(n, 0, 1)

#get roots
roots = roots2(a, b, c)
#booleanize
rootbool = roots>0
rootbool = rootbool[complete.cases(rootbool),]

P = matrix(0, nrow = 2, ncol = 2)
P[1,1] = countevents(rootbool, TRUE, TRUE)
P[2,2] = countevents(rootbool, FALSE, FALSE)
P[1,2] = countevents(rootbool, TRUE, FALSE)
P[2,1] = countevents(rootbool, FALSE, TRUE)
