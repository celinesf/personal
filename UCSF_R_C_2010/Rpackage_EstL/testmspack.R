
R CMD roxygen -d mspack


 package.skeleton('hmstry', code_files='ms.R',force=TRUE) 

package.skeleton(name="EstL", 
code_files=c(
"check_fun.R",
"cmd_line.R",
"estimate.R", 
"get_loci.R",
"get_rho.R",
"ms.R",
"output_cmd.R",
"simulate_data.R", 
 "stats_pop.R", 
"str_fun.R"
)
,force=TRUE) 




 package.skeleton(name="mspack", code_files=c("ms.R","mspack-package.R"),force=TRUE)
 
 library(roxygen) 
roxygenize('mspack',roxygen.dir='mspack',copy.package=FALSE,unlink.target=FALSE) 

library("mspack")
?ms
--------------------
library(Rd2roxygen) 
parse_and_save("stats_pop.Rd",file="stat")
parse_and_save("write_ms_output.Rd",file="write")

-------------------
library("Rmspack")

input_estlik(paramfile="estim2",regfile="regions", datafile="reg0", outdata=TRUE,outparam=TRUE,locifile="loci");


paramfile="estim"
 datafile="reg2"
  regfile="regions"
  locifile="loci"
  #summaries="S_1 S_2 S_sl S_sh F_st D_1 D_2 D_star1 D_star2  D_prime D_prime1 D_prime2, r_square r_square1  	r_square2 nH nH_1 nH_2 Rm Rm_1 Rm_2"
outparam=TRUE
outdata=TRUE
#type=c(0,1,5,10)




---------------
simulate_data(paramfile="sim2", regfile="regions",locifile="loci",datafile="reg2");

simulate_data(paramfile="sim2", regfile="regions",locifile="loci",datafile="reg0",append=FALSE);


paramfile="sim2"
regfile="regions"
datafile="reg2"
#type=c(0,1,5,10)
append=TRUE
locifile="loci"
-----------------------
-------------- order_param
order_param(param, listparam)

-----------
check_param(param=ret$vparam, paramfile=paramfile, listparam=listparam,type=ret$type)

param=ret$vparam
paramfile=paramfile
 listparam=listparam
 type=ret$type


---------------- roxygen
R CMD roxygen 

install.packages('rogygen',type="source")

 library(roxygen) 
 package.skeleton('helloRoxygen', code_files='hello-roxygen.R',force=TRUE) 

R CMD roxygen -d helloRoxygen works, too. 
 roxygenize('helloRoxygen',roxygen.dir='helloRoxygen',copy.package=FALSE,unlink.target=FALSE) 
