param_est=read.csv("param_est.csv",header =FALSE,sep="\n")
param_sim=read.csv("param_sim.csv",header =FALSE,sep="\n")

info_region=read.csv("info_region.csv",header =FALSE)
info_loci=read.csv("info_loci.csv",header =FALSE)
data_file=read.csv("data_file.csv",header =FALSE, blank.lines.skip=FALSE)


list_stats=list()
list_stats[[1]]=c("S","S_1","S_2","S_s","S_ss","S_sf1","S_sf2","S_sl","S_sh","S_f","S_f1","S_f2","S_o","S_o1","S_o2","F(S)","F(S_1)","F(S_2)","F(S_s)","F(S_ss)","F(S_sf1)","F(S_sf2)","F(S_sl)","F(S_sh)","F_st","H_b","H_w1","H_w2","Snn","pi","pi_1","pi_2","theta_W","theta_W1","theta_W2","theta_H","theta_H1","theta_H2","D","D_1","D_2","H","H_1","H_2","D_star","D_star1","D_star2")
#list_stats=read.table("list_stats.csv")
