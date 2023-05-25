#load data

data = read.csv("../data/expected_n_roots_deg_15_div_1_s_1000_normal_stabtrue.csv", sep = "", 
                       header = T)
names(data) = c('d', 'n', 'sim', 'real', 'positive', 'max_eig')

var_curve = function(v_inf, n, d){
  return(v_inf*d^(n/2))
}


df = data %>% group_by(d, n, sim) %>% slice_min(max_eig) 

df_var = df %>% group_by(d, n) %>% 
  summarise(var = sd(positive)^2,
            var_t = mean(positive^2)-(mean(positive))^2) %>% 
  mutate(v_inf = var/sqrt(d),
         var_pred = var_curve(v_inf, n, d)) %>% 
  ungroup() %>% 
  mutate(v_inf_av = mean(v_inf))

ggplot(df_var)+ 
  geom_point(aes(x = d, y = var))+
  geom_point(aes(x = d, y = v_inf),
             color = 'blue',
             shape = '+',
             size = 5)

ggplot(results)+ 
  geom_point(aes(x = d-1, y = obs_variance_real ,
                 color = as.factor(n)))+
  geom_point(aes(x = d-1, y = v_inf),
             color = 'blue',
             shape = '+',
             size = 5)

#plot the variance of my simulations and compare with the predicted variance

var_obs = results$obs_variance_real
var_pred = results$var_pred

plot(var_obs, var_pred)
