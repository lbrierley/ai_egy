# Species level prevalences

mn_data %>%  
  group_by(Species, Result) %>% 
  summarise(n = n()) %>% 
  spread(Result, n) %>%
  replace(is.na(.), 0) %>%
  mutate(prev = positive/(positive+negative)) %>%
  arrange(-prev)

bvbrc_data %>%   
  group_by(Species, Result) %>% 
  summarise(n = n()) %>% 
  spread(Result, n) %>%
  replace(is.na(.), 0) %>%
  mutate(prev = positive/(positive+negative)) %>%
  arrange(-prev)
