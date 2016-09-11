###############################NHANES_d
#####Load Demographic
##############################
Demo <- read.csv ("demo_d.csv", header=T);
Age=cbind(Demo$SEQN,Demo$RIDAGEYR);
Gender=cbind(Demo$SEQN,Demo$RIAGENDR);
Race=cbind(Demo$SEQN,Demo$RIDRETH1);
Data=data.frame(list(SEQN=Demo$SEQN,Wave=rep('d',nrow(Demo))));
Waveid='c';
Demo <- read.csv (paste('demo_',Waveid,'.csv',sep=""), header=T);
Age=rbind(Age,cbind(Demo$SEQN,Demo$RIDAGEYR));
Gender=rbind(Gender,cbind(Demo$SEQN,Demo$RIAGENDR));
Race=rbind(Race,cbind(Demo$SEQN,Demo$RIDRETH1));
Data=rbind(Data,data.frame(list(SEQN=Demo$SEQN,Wave=rep(Waveid,nrow(Demo)))));
Waveid='b';
Demo <- read.csv (paste('demo_',Waveid,'.csv',sep=""), header=T);
Age=rbind(Age,cbind(Demo$SEQN,Demo$RIDAGEYR));
Gender=rbind(Gender,cbind(Demo$SEQN,Demo$RIAGENDR));
Race=rbind(Race,cbind(Demo$SEQN,Demo$RIDRETH1));
Data=rbind(Data,data.frame(list(SEQN=Demo$SEQN,Wave=rep(Waveid,nrow(Demo)))));
Waveid='a';
Demo <- read.csv (paste('demo_',Waveid,'.csv',sep=""), header=T);
Age=rbind(Age,cbind(Demo$SEQN,Demo$RIDAGEYR));
Gender=rbind(Gender,cbind(Demo$SEQN,Demo$RIAGENDR));
Race=rbind(Race,cbind(Demo$SEQN,Demo$RIDRETH1));
Data=rbind(Data,data.frame(list(SEQN=Demo$SEQN,Wave=rep(Waveid,nrow(Demo)))));
#################################
###Load BMI-Waist
#################################
BMX <- (read.csv ("bmx_d.csv", header=T));
BMIC=cbind(BMX$SEQN,as.numeric(as.character(BMX$BMXBMI)));
Waist=cbind(BMX$SEQN,as.numeric(as.character(BMX$BMXWAIST)));
Waveid='c';
BMX <- read.csv (paste('bmx_',Waveid,'.csv',sep=""), header=T);
BMIC=rbind(BMIC,cbind(BMX$SEQN,as.numeric(as.character(BMX$BMXBMI))));
Waist=rbind(Waist,cbind(BMX$SEQN,as.numeric(as.character(BMX$BMXWAIST))));
Waveid='b';
BMX <- read.csv (paste('bmx_',Waveid,'.csv',sep=""), header=T);
BMIC=rbind(BMIC,cbind(BMX$SEQN,as.numeric(as.character(BMX$BMXBMI))));
Waist=rbind(Waist,cbind(BMX$SEQN,as.numeric(as.character(BMX$BMXWAIST))));
Waveid='a';
BMX <- read.csv (paste('bmx_',Waveid,'.csv',sep=""), header=T);
BMIC=rbind(BMIC,cbind(BMX$SEQN,as.numeric(as.character(BMX$BMXBMI))));
Waist=rbind(Waist,cbind(BMX$SEQN,as.numeric(as.character(BMX$BMXWAIST))));
#################################
#################################
#Load Diabetes
#################################
DIQ<-read.csv ("diq_d.csv", header=T);table(DIQ$DIQ010)
T2D=cbind(DIQ$SEQN,DIQ$DIQ010);
Waveid='c';
DIQ <- read.csv (paste('diq_',Waveid,'.csv',sep=""), header=T);
T2D=rbind(T2D,cbind(DIQ$SEQN,DIQ$DIQ010));
Waveid='b';
DIQ <- read.csv (paste('diq_',Waveid,'.csv',sep=""), header=T);
T2D=rbind(T2D,cbind(DIQ$SEQN,as.numeric(as.character(DIQ$DIQ010))));
Waveid='a';
DIQ <- read.csv (paste('diq_',Waveid,'.csv',sep=""), header=T);
T2D=rbind(T2D,cbind(DIQ$SEQN,as.numeric(as.character(DIQ$DIQ010))));
#################################
#Cholestrol
#################################
#Cholestrol:Total,HDL,Total/HDL
#Wave d
TCHOL<-read.csv ("tchol_d.csv", header=T);
TCHOL2<-read.csv ("hdl_d.csv", header=T);
TCH=cbind(TCHOL$SEQN,as.numeric(as.character(TCHOL$LBXTC)),
as.numeric(as.character(TCHOL2$LBDHDD)),
as.numeric(as.character(TCHOL$LBXTC))/as.numeric(as.character(TCHOL2$LBDHDD)));
#Wave c
TCHOL<-read.csv ("l13_c.csv", header=T);
TCH=rbind(TCH,cbind(TCHOL$SEQN,as.numeric(as.character(TCHOL$LBXTC))
,as.numeric(as.character(TCHOL$LBXHDD)),
as.numeric(as.character(TCHOL$LBXTC))/as.numeric(as.character(TCHOL$LBXHDD))));
#Wave b
TCHOL<-read.csv ("l13_b.csv", header=T);
TCH=rbind(TCH,cbind(TCHOL$SEQN,as.numeric(as.character(TCHOL$LBXTC))
,as.numeric(as.character(TCHOL$LBDHDL)),
as.numeric(as.character(TCHOL$LBXTC))/as.numeric(as.character(TCHOL$LBDHDL))));
#Wave a
TCHOL<-read.csv ("lab13_a.csv", header=T);
TCH=rbind(TCH,cbind(TCHOL$SEQN,as.numeric(as.character(TCHOL$LBXTC))
,as.numeric(as.character(TCHOL$LBDHDL)),
as.numeric(as.character(TCHOL$LBXTC))/as.numeric(as.character(TCHOL$LBDHDL))));

