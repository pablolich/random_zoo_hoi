upper_bound_pf <- read.table("../data/p_feas_upper_bound.csv", 
                                       sep = " ", 
                                       header = F)
columnames(upper_bound_pf) = c("n", "d", "pf")
