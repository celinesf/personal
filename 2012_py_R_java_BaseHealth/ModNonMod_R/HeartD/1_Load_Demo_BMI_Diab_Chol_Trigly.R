#################################
# Load CDC Data
#################################
### Load Demographic
Demo <- read.csv ("demo_d.csv", header=T);
Age=cbind(Demo$SEQN,Demo$RIDAGEYR);
Gender=cbind(Demo$SEQN,Demo$RIAGENDR);
Race=cbind(Demo$SEQN,Demo$RIDRETH1);
Data=data.frame(list(SEQN=Demo$SEQN,Wave=rep('d',nrow(Demo))));
Demo <- read.csv ("demo_c.csv", header=T);
Age=rbind(Age,cbind(Demo$SEQN,Demo$RIDAGEYR));
Gender=rbind(Gender,cbind(Demo$SEQN,Demo$RIAGENDR));
Race=rbind(Race,cbind(Demo$SEQN,Demo$RIDRETH1));
Data=rbind(Data,data.frame(list(SEQN=Demo$SEQN,Wave=rep('c',nrow(Demo)))));
Demo <- read.csv ("demo_b.csv", header=T);
Age=rbind(Age,cbind(Demo$SEQN,Demo$RIDAGEYR));
Gender=rbind(Gender,cbind(Demo$SEQN,Demo$RIAGENDR));
Race=rbind(Race,cbind(Demo$SEQN,Demo$RIDRETH1));
Data=rbind(Data,data.frame(list(SEQN=Demo$SEQN,Wave=rep('b',nrow(Demo)))));
Demo <- read.csv ("demo_a.csv", header=T);
Age=rbind(Age,cbind(Demo$SEQN,Demo$RIDAGEYR));
Gender=rbind(Gender,cbind(Demo$SEQN,Demo$RIAGENDR));
Race=rbind(Race,cbind(Demo$SEQN,Demo$RIDRETH1));
Data=rbind(Data,data.frame(list(SEQN=Demo$SEQN,Wave=rep('a',nrow(Demo)))));

### Load BMI-Waist
BMX <- (read.csv ("bmx_d.csv", header=T));
BMIC=cbind(BMX$SEQN,as.numeric(as.character(BMX$BMXBMI)));
Waist=cbind(BMX$SEQN,as.numeric(as.character(BMX$BMXWAIST)));
BMX <- read.csv ("bmx_c.csv", header=T);
BMIC=rbind(BMIC,cbind(BMX$SEQN,as.numeric(as.character(BMX$BMXBMI))));
Waist=rbind(Waist,cbind(BMX$SEQN,as.numeric(as.character(BMX$BMXWAIST))));
BMX <- read.csv ("bmx_b.csv", header=T);
BMIC=rbind(BMIC,cbind(BMX$SEQN,as.numeric(as.character(BMX$BMXBMI))));
Waist=rbind(Waist,cbind(BMX$SEQN,as.numeric(as.character(BMX$BMXWAIST))));
BMX <- read.csv ("bmx_a.csv", header=T);
BMIC=rbind(BMIC,cbind(BMX$SEQN,as.numeric(as.character(BMX$BMXBMI))));
Waist=rbind(Waist,cbind(BMX$SEQN,as.numeric(as.character(BMX$BMXWAIST))));

### Load Diabetes
DIQ<-read.csv ("diq_d.csv", header=T);
T2D=cbind(DIQ$SEQN,DIQ$DIQ010);
DIQ <- read.csv ("diq_c.csv", header=T);
T2D=rbind(T2D,cbind(DIQ$SEQN,DIQ$DIQ010));
DIQ <- read.csv ("diq_b.csv", header=T);
T2D=rbind(T2D,cbind(DIQ$SEQN,as.numeric(as.character(DIQ$DIQ010))));
DIQ <- read.csv ("diq_a.csv", header=T);
T2D=rbind(T2D,cbind(DIQ$SEQN,as.numeric(as.character(DIQ$DIQ010))));

### Load Cholestrol (Total,HDL,Total/HDL)
TCHOL<-read.csv ("tchol_d.csv", header=T);
TCHOL2<-read.csv ("hdl_d.csv", header=T);
TCH=cbind(TCHOL$SEQN,as.numeric(as.character(TCHOL$LBXTC)), as.numeric(as.character(TCHOL2$LBDHDD)), as.numeric(as.character(TCHOL$LBXTC))/as.numeric(as.character(TCHOL2$LBDHDD)));
TCHOL<-read.csv ("l13_c.csv", header=T);
TCH=rbind(TCH,cbind(TCHOL$SEQN,as.numeric(as.character(TCHOL$LBXTC)), as.numeric(as.character(TCHOL$LBXHDD)), as.numeric(as.character(TCHOL$LBXTC))/as.numeric(as.character(TCHOL$LBXHDD))));
TCHOL<-read.csv ("l13_b.csv", header=T);
TCH=rbind(TCH,cbind(TCHOL$SEQN,as.numeric(as.character(TCHOL$LBXTC)), as.numeric(as.character(TCHOL$LBDHDL)), as.numeric(as.character(TCHOL$LBXTC))/as.numeric(as.character(TCHOL$LBDHDL))));
TCHOL<-read.csv ("lab13_a.csv", header=T);
TCH=rbind(TCH,cbind(TCHOL$SEQN,as.numeric(as.character(TCHOL$LBXTC)), as.numeric(as.character(TCHOL$LBDHDL)), as.numeric(as.character(TCHOL$LBXTC))/as.numeric(as.character(TCHOL$LBDHDL))));

### Load Triglycride
TRIGLY<-read.csv ("trigly_d.csv", header=T);
TRIG=cbind(TRIGLY$SEQN,as.numeric(as.character(TRIGLY$LBXTR)));
TRIGLY<-read.csv ("l13am_c.csv", header=T);
TRIG=rbind(TRIG,cbind(TRIGLY$SEQN,as.numeric(as.character(TRIGLY$LBXTR))));
TRIGLY<-read.csv ("l13am_b.csv", header=T);
TRIG=rbind(TRIG,cbind(TRIGLY$SEQN,as.numeric(as.character(TRIGLY$LBXTR))));
TRIGLY<-read.csv ("lab13am_a.csv", header=T);
TRIG=rbind(TRIG,cbind(TRIGLY$SEQN,as.numeric(as.character(TRIGLY$LBXTR))));


#################################
# Processing
#################################
Gender[which(Gender[,2]==1),2]='M';
Gender[which(Gender[,2]==2),2]='F';

Race[which(Race[,2]==1),2]='2Hisp'; # Mexican-American
Race[which(Race[,2]==2),2]='2Hisp'; # Other Hispanic
Race[which(Race[,2]==3),2]='1NHWhite';
Race[which(Race[,2]==4),2]='3NHBlack';
Race[which(Race[,2]==5),2]='4Other';

# BMI as a categorical variable
BMI=BMIC;
BMI[which(BMIC[,2]<18.5),2]='2Underweight';
BMI[which(BMIC[,2]>=18.5&BMIC[,2]<25),2]='1Normal';
BMI[which(BMIC[,2]>=25&BMIC[,2]<30),2]='3Overweight';
BMI[which(BMIC[,2]>=30),2]='4Obese';

TCHL=TCH;
TCHL[which(TCH[,2]<200),2]='1Normal';
TCHL[which(TCH[,2]>=200&TCH[,2]<240),2]='2High';
TCHL[which(TCH[,2]>=240),2]='3VeryHigh';

TRIGL=TRIG;
TRIGL[which(TRIG[,2]<200),2]='1Normal';
TRIGL[which(TRIG[,2]>=200&TRIG[,2]<500),2]='2High';
TRIGL[which(TRIG[,2]>=500),2]='3VeryHigh';

HDL=TCH[,c(1,3)];
HDLL=HDL;
HDLL[which(HDL[,2]>=40),2]='1Normal';
HDLL[which(HDL[,2]<40),2]='2Low';

CHR=TCH[,c(1,4)];
CHRL=HDL;
CHRL[which(CHR[,2]>=9),2]='3VeryHigh';
CHRL[which(CHR[,2]>=6&CHR[,2]<9),2]='2High';
CHRL[which(CHR[,2]<6),2]='1Normal';

T2D[which(T2D[,2]>2),2]='.';
T2D[which(T2D[,2]==2),2]=0;


#################################
# Matching
#################################
loc=match(Data$SEQN,Age[,1]);
Data$Age=as.factor(round(Age[loc,2]/10));
loc=match(Data$SEQN,Age[,1]);
Data$AgeC=as.numeric(Age[loc,2]);
loc=match(Data$SEQN,Gender[,1]);
Data$Gender=as.factor(Gender[loc,2]);
loc=match(Data$SEQN,Race[,1]);
Data$Race=as.factor(Race[loc,2]);
loc=match(Data$SEQN,TCH[,1]);
Data$TCH=as.numeric(TCH[loc,2]);
loc=match(Data$SEQN,HDL[,1]);
Data$HDL=as.numeric(HDL[loc,2]);
loc=match(Data$SEQN,CHR[,1]);
Data$CHR=as.numeric(CHR[loc,2]);
loc=match(Data$SEQN,TCHL[,1]);
Data$TCHL=as.factor(TCHL[loc,2]);
loc=match(Data$SEQN,HDLL[,1]);
Data$HDLL=as.factor(HDLL[loc,2]);
loc=match(Data$SEQN,CHRL[,1]);
Data$CHRL=as.factor(CHRL[loc,2]);
loc=match(Data$SEQN,BMI[,1]);
Data$BMI=as.factor(BMI[loc,2]);
loc=match(Data$SEQN,BMIC[,1]);
Data$BMIC=as.numeric(BMIC[loc,2]);
loc=match(Data$SEQN,Waist[,1]);
Data$Waist=as.numeric(Waist[loc,2]);
loc=match(Data$SEQN,T2D[,1]);
Data$T2D=as.numeric(T2D[loc,2]);
loc=match(Data$SEQN,TRIG[,1]);
Data$TRIG=as.numeric(TRIG[loc,2]);
loc=match(Data$SEQN,TRIGL[,1]);
Data$TRIGL=as.factor(TRIGL[loc,2]);




##########################################
### Modeling for factors corelation
##########################################
WaistModel=lm(Waist~BMIC,data=Data);
Data$WaistAdBMI=(Data$Waist-Data$BMIC*WaistModel$coefficients[2]-WaistModel$coefficients[1]);
CholModel=lm(TCH~BMIC,data=Data);
Data$CholAdBMI=(Data$TCH-Data$BMIC*CholModel$coefficients[2]-CholModel$coefficients[1]);
CholModel=lm(TCH~AgeC,data=Data);


#########################################
### Limiting continous varaibles
#########################################
Choles_ratio_max=13;
Choles_ratio_min=5;
BMI_max=35;
BMI_min=22.5;
Waist_min=-20;
Waist_max=20;
TRIG_max=500;
TRIG_min=100;

Data$CHR2=Data$CHR;
Data[which(Data$CHR2>Choles_ratio_max),]$CHR2=Choles_ratio_max;
Data[which(Data$CHR2<Choles_ratio_min),]$CHR2=Choles_ratio_min;

Data$BMIC2=Data$BMIC;
Data[which(Data$BMIC2>BMI_max),]$BMIC2=BMI_max;
Data[which(Data$BMIC2<BMI_min),]$BMIC2=BMI_min;

Data$WaistAdBMI2=Data$WaistAdBMI;
Data[which(Data$WaistAdBMI2>Waist_max),]$WaistAdBMI2=Waist_max;
Data[which(Data$WaistAdBMI2<Waist_min),]$WaistAdBMI2=Waist_min;



Data$TRIG2=Data$TRIG;
Data[which(Data$TRIG2>TRIG_max),]$TRIG2=TRIG_max;
Data[which(Data$TRIG2<TRIG_min),]$TRIG2=TRIG_min;


Data$CholAdAge=(Data$TCH-Data$AgeC*CholModel$coefficients[2]-CholModel$coefficients[1]);