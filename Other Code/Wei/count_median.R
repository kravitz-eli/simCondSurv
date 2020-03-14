# See how many posterior predictive survival times reached 50% survival -------
time_and_event[[1]] %>% 
  group_by(trt) %>% 
  summarise(n_event = sum(event))

event_count = map2_dfr(
  time_and_event, 
  1:length(time_and_event),
  ~ {.x %>% 
      group_by(trt) %>% 
      summarise(n_event = sum(event)) %>% 
      mutate("posterior_num" = .y) %>% 
      select(posterior_num, everything())
  })

event_count %>% 
  mutate("median_exists" = n_event > 334 / 2) %>% 
  group_by(trt) %>% 
  summarise("count_median_exists" = sum(median_exists))

# See median survival time in trt arm when we observe 300 events
cens_and_median = time_and_event %>% 
  map_dfr(., ~{
    .x %>%
      filter(trt == 1) %>% 
      summarise("median_surv" = median(time),
                "n_censor" = sum(1 - event)) 
  })
