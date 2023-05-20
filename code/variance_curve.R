#load data

data = read.csv("../data/expected_n_roots_dim_25_div_25_s_10000_normal_stabtrue.csv", sep = "", 
                       header = T)
names(data) = c('d', 'n', 'sim', 'real', 'positive', 'max_eig')

df = data %>% group_by(d, n, sim) %>% slice_min(max_eig) 

df_var = df %>% group_by(d, n) %>% summarise(var = sd(positive)^2) %>% 
  mutate(v_inf = 4*var/sqrt(d)) %>% 
  ungroup() %>% 
  mutate(v_inf_av = mean(v_inf))

ggplot(df_var)+ 
  geom_point(aes(x = d, y = var))+
  geom_point(aes(x = d, y = v_inf),
             color = 'blue')

#plot the variance of my simulations and compare with the predicted variance

var_obs = results$obs_variance_pos
var_pred = results$var_pred

plot(var_obs, var_pred)
