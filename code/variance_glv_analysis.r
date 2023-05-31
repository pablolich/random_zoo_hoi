require(tidyverse)
vars = read.csv("../data/variances.csv", header = F, sep="\t")
h = hist(vars$V1, breaks = seq(max(vars$V1)),
         main = "monomial var n=2 d=4")

variance_df=tibble()
for (h in seq(2, 10)){
  for (p in seq(1,9)){
    if (p>=h){
      variance_df = rbind(variance_df, tibble(
        p = p,
        h = h,
        var = -1)
      )
    }else{
      variance_df = rbind(variance_df, tibble(
        p = p,
        h = h,
        var = factorial(h-1)/(factorial(p)*factorial(h-1-p))
      ))
    }
  }
}

#save data
write.table(variance_df, '../data/variance_structure.dat', sep = " ", row.names = F,
            quote = F)
