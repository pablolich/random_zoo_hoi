require(tidyverse)
vars = read.csv("../data/variances.csv", header = F, sep="\t")
h = hist(vars$V1, breaks = seq(max(vars$V1)),
         main = "monomial var n=2 d=4")
