suppressMessages({library(macrosheds); library(dplyr)})
fc <- macrosheds::file_ids_for_r_package
print(names(fc))
doms <- c("fernow", "hbef", "hjandrews", "santa_barbara", "santee")
print(as.data.frame(fc %>% filter(domain %in% doms)))
