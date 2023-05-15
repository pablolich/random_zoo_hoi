library("tidyverse")
library("viridis")

expected_eq = function(n, d = 3){
  return (1/(2^(n))*(d-1)^((n)/2))
}

#when load new data, delete first row of columns!

data = read.csv("../data/expected_n_rootsdim_6_div_6_s_2000.csv", sep = "", 
                 header = F )

names(data) = c('d', 'n', 'sim', 'n_positive')

n_sim = 2000
  
df = data %>% group_by(d, n) %>% 
  count(positive, name = 'count') %>% 
  mutate(total_eq = sum(positive*count),
         mean_eq = total_eq/n_sim,
         var_eq = mean((positive-mean_eq)^2))

  ggplot(df, aes(x = n, y = exp_eq))+
    geom_point(aes(color = as.factor(d)))+
    geom_function(fun = expected_eq, args = list(d = 2),
                  lty = 2)+
    geom_function(fun = expected_eq, args = list(d = 3))+
    geom_function(fun = expected_eq, args = list(d = 4))+
    geom_function(fun = expected_eq, args = list(d = 5))+
    geom_function(fun = expected_eq, args = list(d = 6))+
    theme(aspect.ratio = 1)

#check feasibility hypothesis
  
data_feas = read.csv("../data/expected_n_rootsdim_4_div_4_s_500.csv",
                     sep = "", 
                     header = F)
names(data_feas) = c('d', 'n', 'sim', 'real', 'positive')
    
n_sim = 500
df = data_feas %>% group_by(d, n) %>% 
  count(real, positive, name = 'n_roots') %>% 
  filter(!(real == 0 & positive ==0)) %>% 
  mutate(tot_real = sum(real*n_roots),
         tot_positive = sum(positive*n_roots))

ggplot(df, aes(x = n, y = tot_positive))+
  geom_point(aes(color = as.factor(d)))+
  geom_point(aes(x = n, y = tot_real/2^n),
             shape = 3)+
  theme(aspect.ratio = 1)

#stability hypothesis

data_stab = read.csv("../data/expected_n_roots_dim_3_div_3_s_10_normal.csv",
                     sep = "", 
                     header = F)
       
             