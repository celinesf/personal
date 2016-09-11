paramfile=read.csv("paramfile.csv",header =FALSE,sep="\n")
listparam=c("theta_1","M_present","theta_2","theta_A","T_split","T_change","M_change","rho", "nregions","type","stats","howmany", "parallel" )

regfile=read.csv("regfile.csv",header =FALSE)
#locifile=read.csv("locifile.csv",header =FALSE)
msout=read.csv("msout.csv",header =FALSE, blank.lines.skip=FALSE)


liststats=list()
liststats[[1]]=c("S","S_1","S_2","S_s","S_ss","S_sf1","S_sf2","S_sl","S_sh","S_f","S_f1","S_f2","S_o","S_o1","S_o2","F(S)","F(S_1)","F(S_2)","F(S_s)","F(S_ss)","F(S_sf1)","F(S_sf2)","F(S_sl)","F(S_sh)","F_st","H_b","H_w1","H_w2","Snn","pi","pi_1","pi_2","theta_W","theta_W1","theta_W2","theta_H","theta_H1","theta_H2","D","D_1","D_2","H","H_1","H_2","D_star","D_star1","D_star2")

