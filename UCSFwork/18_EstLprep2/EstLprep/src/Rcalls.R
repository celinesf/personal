
library("Rmspack")

m=get_ms_output("m1")
stats_pop(m,type = c(0, 1, 5, 10))

stats_pop(m,type = c(1, 1, 5, 10))