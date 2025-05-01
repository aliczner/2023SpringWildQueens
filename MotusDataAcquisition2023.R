##getting motus detections

library(motus)
library(dplyr)
library(lubridate)

#function tagme downloads data from project number x (437)
sql.motus <- tagme(437, new = TRUE, update = TRUE)

tbl.alltags <- tbl(sql.motus, "alltags") #retrieve alltags from the SQLite database

df.alltags <- tbl.alltags %>% #convert to dataframe
  collect() %>%
  mutate(time = as_datetime(ts))

df.selected <- df.alltags %>%
  select(22, 33, 45, 46, 50, 63)


df.2023 <- df.selected %>%
  filter(year(time) == 2023)

df.spring <- df.2023 %>%
  filter(month(time) %in% c(4, 5, 6)) #6316 rows

write.csv(df.spring, "2023SpringMotus.csv")

df.summer <- df.2023 %>%
  filter(month(time) %in% c(8, 9, 10, 11)) #34969 rows

write.csv(df.summer, "2023SummerMotus.csv")
