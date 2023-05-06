library("tidyverse")
library("viridis")

expected_eq = function(n, d = 3){
  return (1/(2^(n))*(d-1)^((n)/2))
}

data = read.csv("../data/expected_n_roots.csv", sep = "", header = F )

n_sim = 1000
  
df = data %>% group_by(V1) %>% 
  count(V3) %>% 
  summarise(exp_eq = sum(V3*n)/n_sim)

  ggplot(df, aes(x = V1, y = exp_eq))+
    geom_point()+
    geom_function(fun = expected_eq, args = list(d = 3))
  
