library("tidyverse")
library("viridis")

expected_eq = function(n, d){
  return ((d)^(n/2)/(2^(n)))
}

f_x = function(n, i){
  return((1-1/2^n)^i)
}

f2_x = function(n, i){
  return(log(1-1/2^n)^2*(1-1/2^n)^(i))
}

f3_x = function(n, i){
  return((1 - 2^(-n))^i*log(1 - 2^(-n))^3)
}

expected_value_expansion = function(n, d, var, skw, order){
  mean_x = sqrt(d^n)
  #evaluate function and second derivative at mean
  f_mean = f_x(n, mean_x)
  f2_mean = f2_x(n, mean_x)
  f3_mean = f2_x(n, mean_x)
  if (order == 1){
    approx = f_mean 
  } else if (order == 2){
    approx = f_mean + 1/2*f2_mean*var 
  } else{
    approx = f_mean + 1/2*f2_mean*var + 1/6*f3_mean*skw
  }
  return(approx)
} 

data = read.csv("../data/expected_n_rootsdim_6_div_6_s_2000.csv", sep = "", 
                 header = F )

names(data) = c('d', 'n', 'sim', 'positive')

n_sim = 2000
  
df = data %>% 
  group_by(d, n) %>%
  count(positive, name = 'count') %>% 
  mutate(total_eq = sum(positive*count),
         mean_eq = total_eq/n_sim,
         var_eq = sum(count*(positive - mean_eq)^2)/n_sim,
         pfeas = 1 - expected_value_expansion(n, d, NA, NA, 1)) %>% 
  arrange(d, n)

  ggplot(df, aes(x = n, y = mean_eq))+
    geom_point(aes(color = as.factor(d)))+
    # geom_errorbar(aes(x = n, y = mean_eq, 
    #                   ymin = mean_eq - sqrt(var_eq),
    #                   ymax = mean_eq + sqrt(var_eq),
    #                   color = as.factor(d)),
    #               size = 0.5,
    #               width = 0.2)+
    geom_point(aes(x = n, y = expected_eq(n, d-1)),
               shape=3)+
    theme(aspect.ratio = 1)

#check feasibility hypothesis

df_feas = data %>% group_by(d, n) %>% 
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

###timing

data = read.csv("../data/expected_n_roots_deg_5_div_5_s_1_normal_stabfalse.csv", sep = "", 
                header = T )

names(data) = c('d', 'n', 'sim', 'real', 'positive', 'time')

data$exp_time = data$time*5000
data$exp_time_h = data$exp_time/3600

#compare to stefano's

data_sa = results %>% select(d, n, time_sec)

merge(data, data_sa, by=c("d", "n"))
