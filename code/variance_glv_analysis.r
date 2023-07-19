require(tidyverse)
library(gtools)
vars = read.csv("../data/variances.csv", header = F, sep="\t")
h = hist(vars$V1, breaks = seq(max(vars$V1)),
         main = "monomial var n=2 d=4")

is_row_there = function(row, df){
  for (i in seq(nrow(df))){
    row_i = df[i,]
    are_equal= prod(row_i == row)
    if (are_equal){
      return(T)
    }
  }
  return(F)
}

get_unique = function(df){
  #get length
  df_unique = tibble()
  number_rows = nrow(df)
  for (i in seq(number_rows)){
    #get row i
    df_i = df[i,]
    unique_row = sort(df_i)
    if (i==1){
      df_unique = rbind(df_unique, unique_row)
    }
    #check if I have seen it in df_unique
    is_present = is_row_there(unique_row, df_unique)
    if (!is_present){
      df_unique = rbind(df_unique, unique_row)
    }
  }
  return(df_unique)
}


variance_as_kss = function(d, l_vec){
  return(factorial(d)/(prod(c(factorial(l_vec), factorial(d-sum(l_vec))))))
}

variance_as_symmetric = function(l_vec){
  return((factorial(sum(l_vec))/(prod(factorial(l_vec))))^2)
}

degree = 4
  
df = permutations(n=degree, r=degree, v = seq(0, 3), repeats.allowed = T)

df_unique = get_unique(df)

df_l_vecs = df_unique %>%
  group_by(X0L, X0L.1, X0L.2, X0L.3) %>%
  mutate(sum = sum(c(X0L, X0L.1, X0L.2, X0L.3))) %>% 
  filter(sum<=degree) %>% 
  arrange(sum) %>% 
  ungroup()

colnames(df_sum) = c("x1", "x2", "x3", "x4", "sum")



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

variance_a_kss=tibble()
variance_a_symmetric=tibble()
nrows_l_vec = nrow(df_l_vecs)
for (h in seq(0, degree)){
  for (l_i in seq(nrows_l_vec)){
    p = df_l_vecs[l_i,]$sum
    lvec_i = df_l_vecs[l_i,][1:3]
    if (p>h){
      variance_a_kss = rbind(variance_a_kss, tibble(
        p = l_i,
        h = h,
        var = -1)
      )
      variance_a_symmetric = rbind(variance_a_symmetric, tibble(
        p = l_i,
        h = h,
        var = -1)
      )
    }else{
      variance_a_kss = rbind(variance_a_kss, tibble(
        p = l_i,
        h = h,
        var = variance_as_kss(h, as.numeric(lvec_i))
      ))
      variance_a_symmetric = rbind(variance_a_symmetric, tibble(
        p = l_i,
        h = h,
        var = variance_as_symmetric(as.numeric(lvec_i))
      ))
    }
  }
}

#save data
write.table(variance_df, '../data/variance_structure.dat', sep = " ", row.names = F,
            quote = F)

write.table(variance_a_kss, '../data/variance_a_kss.dat', sep = " ", row.names = F,
            quote = F)

write.table(variance_a_symmetric, '../data/variance_a_symmetric.dat', sep = " ", row.names = F,
            quote = F)


