library(tidyverse)


prof_times <- read.csv("./profiling_times.csv")
prof_times_NO <- read.csv("./profiling_times_NO_scs.csv")


prof_times <- prof_times %>% mutate(type = "normal")
prof_times_NO <- prof_times_NO %>% mutate(type = "no_scs")

colums = c("samples", "size", "comp_time", "type_sim")

names(prof_times) <- colums
names(prof_times_NO) <- colums

df <- rbind(prof_times, prof_times_NO)

df$type <- as.factor(df$type_sim)
#df$type <- factor("type", levels = "normal", "no_scs")
df$size <- as.numeric(df$size)

str(df)

ggplot(df, aes(x=size, y=comp_time, color=type_sim)) +
  geom_point(show.legend = F) + geom_smooth(show.legend = F)+
  ylim(c(0, 100000)) +
  scale_x_discrete(limit = c(2, 4, 8, 16, 32)) +
  xlab("Community Size") +
  ylab("Computation Time (s)")+
  theme_bw() +
  geom_label(label = "SCS in Simulation",
             x = 25, y = 75000, color = "black") +
  geom_label(label = "No SCS in Simulation",
             x = 25, y = 10000, color = "black")
  
