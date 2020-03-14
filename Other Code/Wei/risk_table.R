mona_L <- survSplit(Surv(time, event) ~ trt, 
                    data = data,
                    id = "id",
                    cut = seq(2, 36, by = 2),
                    start = "start",
                    episode ="interval") %>% 
  mutate(
    interval = factor(
      interval, 
      labels = paste0("(", seq(0, 34, by = 2), ",", seq(2, 36, by = 2), ")")
    )
  )
as_tibble()


risk_tab = mona_L %>% 
  group_by(interval, trt) %>% 
  summarise("num_at_risk" = sum((1 - event))) %>% 
  pivot_wider(id_cols = interval, 
              values_from = num_at_risk, 
              names_from = trt) %>% 
  select(
    interval,
    "riboc+let" = `1`,
    "plac+let" = `0`
  )


fit <- survfit(Surv(time, event) ~ trt, data = data)
library(survminer)
ggrisktable(
  fit = fit,
  data = data, 
  break.time.by = 2,
  risk.table.type = c("absolute"),
  
  )
