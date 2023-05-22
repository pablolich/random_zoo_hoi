library("tidyverse")
library("viridis")

expected_eq = function(n, d){
  return (1/(2^(n))*(d)^(n/2))
}

#when load new data, delete first row of columns!

  data = read.csv("../data/expected_n_roots_dim_5_div_5_s_1000_normal_stabtrue.csv", sep = "", 
                 header = F )

names(data) = c('d', 'n', 'sim', 'n_positive')

n_sim = 1000
  
df = data %>% group_by(d, n) %>% 
  slice_min(max_eig) %>% 
  count(n_positive, name = 'count') %>% 
  mutate(total_eq = sum(n_positive*count),
         mean_eq = total_eq/n_sim,
         var_eq = sum(count*(n_positive - mean_eq)^2)/n_sim)

  ggplot(df, aes(x = n, y = mean_eq))+
    geom_point(aes(color = as.factor(d)))+
    geom_errorbar(aes(x = n, y = mean_eq, 
                      ymin = mean_eq - sqrt(var_eq),
                      ymax = mean_eq + sqrt(var_eq),
                      color = as.factor(d)),
                  size = 0.5,
                  width = 0.2)+
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

names(data_stab) = c('d', 'n', 'sim', 'real', 'positive', 'max_eig')

#check results with uniform distribution

data_unif = read.csv("../data/expected_n_roots_dim_5_div_5_s_1000_uniform_stabtrue.csv", sep = "", 
                header = T )

data_normal = read.csv("../data/expected_n_roots_dim_5_div_5_s_1000_normal_stabtrue.csv", sep = "", 
                       header = T)
names(data_unif) = c('d', 'n', 'sim', 'real', 'positive', 'max_eig')
names(data_normal) = c('d', 'n', 'sim', 'real', 'positive', 'max_eig')

n_sim = 1000

df = data_normal %>% 
  group_by(d, n, sim) %>% 
  slice_min(max_eig) %>% 
  group_by(d, n) %>% 
  count(positive, name = 'count') %>% 
  mutate(mean_eq = sum(positive*count)/n_sim,
         var_eq = sum(count*(positive - mean_eq)^2)/n_sim)

df_u = data_unif %>% 
  group_by(d, n, sim) %>% 
  slice_min(max_eig) %>% 
  group_by(d, n) %>% 
  count(positive, name = 'count') %>% 
  mutate(mean_eq = sum(positive*count)/n_sim,
         var_eq = sum(count*(positive - mean_eq)^2)/n_sim)

ggplot(df, aes(x = n, y = mean_eq))+
  geom_point(aes(color = as.factor(d)))+
  geom_point(aes(x = n, y = expected_eq(n, d)),
             shape=3)+
  theme(aspect.ratio = 1)+
  labs(title = 'Gaussian')

ggplot(df_u, aes(x = n, y = mean_eq))+
  geom_point(aes(color = as.factor(d)))+
  geom_point(aes(x = n, y = expected_eq(n, d)),
             shape=3)+
  theme(aspect.ratio = 1)+
  labs(title = 'Uniform')