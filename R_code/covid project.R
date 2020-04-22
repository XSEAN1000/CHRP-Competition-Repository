Data = time_series_covid19_confirmed_US
globaldata = time_series_covid19_confirmed_global

####COUNTRIES BESIDES US##################################################
####COUNTRIES BESIDES US##################################################
####COUNTRIES BESIDES US##################################################
####COUNTRIES BESIDES US##################################################
####COUNTRIES BESIDES US##################################################
####COUNTRIES BESIDES US##################################################
####COUNTRIES BESIDES US##################################################

Canada = globaldata[36:46,]
L = ncol(Canada)
Canadats = Canada[,5:L]

Italyts = globaldata[138,5:L]
Israelts = globaldata[137,5:L]
Mexicots = globaldata[159,5:L]
Japants = globaldata[140,5:L]
NZts = globaldata[171,5:L]
Pakistants = globaldata[178,5:L]
Vietnamts = globaldata[229,5:L]
SAfricats = globaldata[201,5:L]

###USA STATES#############################################################
###USA STATES#############################################################
###USA STATES#############################################################
###USA STATES#############################################################
###USA STATES#############################################################
###USA STATES#############################################################
###USA STATES#############################################################
###USA STATES#############################################################
###USA STATES#############################################################
###USA STATES#############################################################


AL = Data[6:72,]
L = ncol(AL)
ALnum = AL[,9:L]
l = ncol(ALnum)
AL_coords = ALnum[,1:2]
ALts = ALnum[,4:l]


CA = Data[192:249,] #all info on CA
CAnum = CA[,9:L] 
CA_coords = CAnum[,1:2]
CAts = CAnum[,4:l]


CO = Data[250:313,] #all info on COloroado
COnum = CO[,9:L] 
CO_coords = COnum[,1:2]
COts = COnum[,4:l]

DC = Data[325,] #all 
DCnum = DC[,9:L]
DC_coords = DCnum[,1:2]
DCts = DCnum[,4:l]

FL = Data[326:392,] #all info on FL
FLnum = FL[,9:L] 
FL_coords = FLnum[,1:2]
FLts = FLnum[,4:l]

GA = Data[393:551,] 
GA_coords = GA[,9:10]
GAts = GA[,12:L]

HI = Data[552:556,] 
HI_coords = HI[,9:10]
HIts = HI[,12:L]

ID = Data[557:600,] 
ID_coords = ID[,9:10]
IDts = ID[,12:L]

IL = Data[601:702,] 
IL_coords = IL[,9:10]
ILts = IL[,12:L]

IN = Data[703:794,] 
IN_coords = IN[,9:10]
INts = IN[,12:L]

IA = Data[795:893,] 
IA_coords = IA[,9:10]
IAts = IA[,12:L]

KS = Data[894:998,] 
KS_coords = KS[,9:10]
KSts = KS[,12:L]

NJ = Data[1780:1800,] 
NJ_coords = NJ[,9:10]
NJts = NJ[,12:L]

NM = Data[1801:1833,] 
NM_coords = NM[,9:10]
NMts = NM[,12:L]

NY = Data[1834:1895,]    #all info on NY
NYnum = NY[,9:L]        #mostly numbers (no names etc)
NY_coords = NYnum[,1:2]  #lat and long
NYts = NYnum[,4:l]      #matrix of TS data


TN = Data[2434:2528,] 
TN_coords = TN[,9:10]
TNts = TN[,12:L]

TX = Data[2529:2782,] 
TX_coords = TX[,9:10]
TXts = TX[,12:L]

plot(TX_coords[,2], TX_coords[,1])

##########################################################################
######WEATHER DATA FROM JAN 1 - current################################
######################################################################

#V2- V5 max, min, avg, departure (temp)

temporary = list.files(pattern = ".txt")
?list.files
weather = lapply(temporary, read.delim, header=FALSE)
names(weather) = c("Ada", "Bern", "Suffolk",
                   "Fairfield", "Clark", "Cook",
                   "Dallas","DC","Denver","Wayne",
                   "Erie", "Fulton", "Harris",
            "LA","Broward", "Maricopa","Miami", "Montgomery",
            "Middle","NY","Jefferson","Philly","Riverside",
            "Sandiego", "Seattle","Sioux", "Wyandotte")

  
LosAngelests = CAts[19,]
NewYorkts = Data[1864,12:L]
Fairfieldts = Data[314,12:L]
Miamits = Data[368,12:L]
Cookts = Data[616,12:L]
Orleansts = Data[1154,12:L]
Jeffersonts = Data[1144,12:L]
Suffolkts = Data[1235,12:L]
Waynets = Data[1318,12:L]
Middlesexts = Data[1791,12:L]
Philts = Data[2300,12:L]
Harrists = Data[2629,12:L]
Kingts = Data[2975,12:L]

Maricopats = Data[109,12:L]
Dallasts = Data[2585,12:L]
Montgomeryts = Data[2698,12:L]
DCts = Data[325,12:L]
Wyandottets = Data[998, 12:L]
Denverts = Data[266,12:L]
Fultonts = Data[452,12:L]
Clarkts = Data[1754, 12:L]
Bernts = Data[1801, 12:L]
Sandiegots = Data[228, 12:L]
Eriets = Data[1848, 12:L]
Adats = Data[557, 12:L]

Browardts = Data[331, 12:L]
Riversidets = Data[224,12:L]
Siouxts = Data[2416,12:L]


###TOOLS#######################################################################
###TOOLS#######################################################################
###TOOLS#######################################################################
###TOOLS#######################################################################
###TOOLS#######################################################################
###TOOLS#######################################################################
###TOOLS#######################################################################
###TOOLS#######################################################################

total <- function(data){      #get case totals for a state, given county data
  N = ncol(data)
  out = matrix(0,nrow = N,ncol = 2)
  for(i in 1:N){
    out[i,2] = sum(data[,i])
    out[i,1] = i
  }
  return(out)
}

shiftm <- function(total,k){ #shifts 'total' datamatrix, for PWE plots
  N = nrow(total)
  temp = total
  for(i in 1:(N%%k)){
    temp = rbind(c(0,0),temp)
  }
  N = nrow(temp)
  for(i in 1:N){
    temp[i,1] = i
  }
  return(temp)
}


expf <- function(data){
  vals = as.vector(total(data)[,2])
  N = length(vals)
  factors = c()
  for(i in 1:(N-1)){
    if(vals[i] == 0){factors[i] = 0}
    else{factors[i] = vals[i+1]/vals[i]}
  }
  return(factors)
}

part <- function(data,k){        #adds zeros to begining of matrix, ncol = k
  N = length(data)
  temp = data
  for(i in 1:(N%%k)){
    temp = c(c(0),temp)
  }
  R = length(temp)/k
  partition = matrix(0,nrow = R,ncol = k)
  for(i in 1:R){
    partition[i,] = temp[1:k]
    temp = temp[-c(1:k)]
  }
  return(partition)
}

part2 <- function(data,k){ #adds zeros at beginning for divisiblity by k
  N = length(data)
  temp = data
  for(i in 1:(N%%k)){
    temp = c(c(0),temp)
  }
  return(temp)
}

nz <- function(data){          #returns nonzero part of vector 
  return(data[data>0])
}

exp_av <- function(data, k){   #gives local approximation by exponentials
  points = total(data)[,2]
  points = part2(points, k)  #added 0 to beginning to make /k
  Q = length(points)/k
  output = matrix(0,nrow = Q, ncol = 2)
  for(i in 1:Q){
    starting = points[1:k]
    if(starting[1]==0){
      output[i,2] = 0
      output[i,1] =0}
    else{
      X = c(1:k) + (i-1)*k
      fit = lm(log(starting) ~ X)
      output[i,2] = fit$coefficients[2]
      output[i,1] = fit$coefficients[1]
    }
    points = points[-c(1:k)]  
  }
  return(output)
}

averager <- function(data, k){   #averages temp to match exponential partition..
  points = data
  points = part2(points, k)
  R = length(points)/k
  output = c() 
  for(i in 1:R){
    starting = points[1:k]
    output[i] = mean(starting)
    points = points[-c(1:k)]  
  }
  return(output)
}

maxer <- function(data, k){   #maxes temp to match exponential partition..
  points = data
  points = part2(points, k)
  R = length(points)/k
  output = c() 
  for(i in 1:R){
    starting = points[1:k]
    output[i] = max(starting)
    points = points[-c(1:k)]  
  }
  return(output)
}

#want to throw away the points before case total = 1
fnz <- function(data){
  L = length(data)
  for(i in 1:L){
    if(data[i] != 0){ return(i)
      break}
  }
}

fn1 <- function(data){     #for getting multiplicative rate first time >0
  L = length(data)
  for(i in 1:L){
    if(data[i] > 1){ return(i)
      break}
  }
}

deriv = function(data, k){  #input cbind(column x, column y)
  if(k == 0){
    x = data[,1] 
    y = data[,2]
  }else{
    n = nrow(data)
    index = c()
    j=1
    while(k*j < n){
      index[j] = k*j
      j = j + 1}
    x = data[index,1]
    y = data[index,2]}
  yprime = c()
  N = length(x)-1
  for(i in 1:N){
    yprime[i] = (y[i+1] - y[i])/ (x[i+1]- x[i])
  }
  return(cbind(x,yprime))
}

dim(weather$LA)

ofset = function(matrix, k){          #offsets temp. data
  a = 22-k
  b = 107- k
  return(matrix[,a:b])
}

BH = function(pvals, q){        #Benjamini-Hochberg 
  N = length(pvals)
  vee = rbind(pvals, c(1:N))
  ordered = vee[,order(vee[1,])]
  indices = c()
  for(i in 1:N){
    if( ordered[1,i] < (q*i)/N){
      indices = c(indices, i)}}
  indices = c(indices,0)
  imax = max(indices)
  if(imax == 0){
    return(c(0))}
  else{return(ordered[2,c(1:imax)])}
}
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################




#####################################################################
#GROWTH RATE APPROXIMATIONS:#########################################
#####################################################################

##PWE:

#Below we plot PWE approximations to case totals, NY or Chicago

temp = NYts
temp = Cookts

par(mfrow = c(2,2))
for(k in 2:5){
  M = exp_av(temp, k)  #says which exp fxn, on 1st k, 2nd k, ...
  L = nrow(M)
  plot(shiftm(total(temp),k), main = "PWE of case total, NY",
      xlab = "day", ylab = "#cases")
  for(i in 1:L){
    X = c(0:k) + (i-1)*k
    lines(X,exp((M[i,2])*X +M[i,1]), col = "red")}
}

#to see what PWE rate/coefficient means, look at:

par(mfrow = c(2,2))
for(k in 2:5){                    #plot c in be^{ct} vs t
  M = exp_av(temp, k)[,2]
  L = length(M)
  plot(1:L,M, main = "PWE. rate coefficient over time")
}

for(k in 2:5){                      #plot b in be^{ct} vs t
  M = exp_av(temp, k)[,1]
  L = length(M)
  plot(1:L,M, main = "PWE. scale coefficient over time")
}

#LINEAR DIFFERENCE


#Below, we demonstrate our "thinning" first difference formula

par(mfrow = c(1,1))
x = seq(0,8,.1)
y = sin(x)
datamatrix = cbind(x,y)

par(mfrow = c(1,1))
plot(datamatrix)
points(deriv(datamatrix,0), col = "red")   #no thinning
points(deriv(datamatrix,2), col = "green4")   #thinning factor 2
points(deriv(datamatrix,5), col = "blue")     #thinning factor 5


#here is the 1st linear difference for case total, Chicago
par(mfrow = c(2,2))
for(k in 2:5){                    #Plot N'(t) vs. t
  M = deriv(total(temp), k)
  plot(M, main = "1st deriv over time")
}


##############################################################################
###########################################################################
#CORRELATION TESTING##################################################
###########################################################################


Places = rbind(LosAngelests,NewYorkts,Fairfieldts,Miamits,Cookts,
               Jeffersonts,Suffolkts,Waynets,Middlesexts,Philts,
               Harrists,Kingts, Maricopats, Dallasts, Montgomeryts,
               DCts, Wyandottets,Denverts,Fultonts,Clarkts,
               Bernts, Sandiegots,Eriets,Adats, Browardts, Riversidets,
               Siouxts
               )
dim(Places)              #27 locations!

TempMaxes = rbind(weather$LA[,2], weather$NY[,2], weather$Fairfield[,2],
                  weather$Miami[,2], weather$Cook[,2], weather$Jefferson[,2],
                  weather$Suffolk[,2],weather$Wayne[,2], weather$Middle[,2],
                  weather$Philly[,2],weather$Harris[,2],weather$Seattle[,2], 
                  weather$Maricopa[,2], weather$Dallas[,2], weather$Montgomery[,2],
                  weather$DC[,2], weather$Wyandotte[,2],weather$Denver[,2],
                  weather$Fulton[,2],weather$Clark[,2],weather$Bern[,2],
                  weather$Sandiego[,2],weather$Erie[,2],weather$Ada[,2],
                  weather$Broward[,2], weather$Riverside[,2],weather$Sioux[,2])

TempMins = rbind(weather$LA[,3], weather$NY[,3], weather$Fairfield[,3],
                 weather$Miami[,3], weather$Cook[,3], weather$Jefferson[,3],
                 weather$Suffolk[,3],weather$Wayne[,3], weather$Middle[,3],
                 weather$Philly[,3],weather$Harris[,3],weather$Seattle[,3], 
                 weather$Maricopa[,3], weather$Dallas[,3], weather$Montgomery[,3],
                 weather$DC[,3], weather$Wyandotte[,3],weather$Denver[,3],
                 weather$Fulton[,3],weather$Clark[,3],weather$Bern[,3],
                 weather$Sandiego[,3],weather$Erie[,3],weather$Ada[,3],
                 weather$Broward[,3],weather$Riverside[,3],weather$Sioux[,3])

TempAvgs = rbind(weather$LA[,4], weather$NY[,4], weather$Fairfield[,4],
                 weather$Miami[,4], weather$Cook[,4], weather$Jefferson[,4],
                 weather$Suffolk[,4],weather$Wayne[,4], weather$Middle[,4],
                 weather$Philly[,4],weather$Harris[,4],weather$Seattle[,4], 
                 weather$Maricopa[,4], weather$Dallas[,4], weather$Montgomery[,4],
                 weather$DC[,4], weather$Wyandotte[,4],weather$Denver[,4],
                 weather$Fulton[,4],weather$Clark[,4],weather$Bern[,4],
                 weather$Sandiego[,4],weather$Erie[,4],weather$Ada[,4],
                 weather$Broward[,4],weather$Riverside[,4],weather$Sioux[,4])

Temps = list(TempMaxes, TempMins, TempAvgs)
Tnames = c("Max T.", "Min T.", "Avg T.")
whichmethod = c("spearman", "kendall", "pearson")
colors = c("red","deeppink","green4", "orange","blue","purple",
           "gold","magenta","green","azure",
           "coral", "cyan", "magenta1","magenta2","magenta3","magenta4",
           "grey",
           "grey","grey","grey","grey","grey","grey",
           "grey","grey","grey","grey")

###################################################################################
###########################################################################
#########################################################################
########MULTIPLICATIVE FACTORS#######################################
#########################################################################
length(expf(Places[1,]))
MF = matrix(0, nrow = 27,ncol = 85)
for(i in 1:27){
  MF[i,] = expf(Places[i,])            #MF data
}
dim(MF)


par(mfrow = c(2,3),mar = c(4,4,4,4))
Z = matrix(0,nrow = 3, ncol = 15)
rhos = list(Z,Z,Z)
pvals = list(Z,Z,Z)
Intercepts = list(Z,Z,Z)
Slopes = list(Z,Z,Z)
x = seq(-5,110,1)
for(l in 1:3){   #which test
for(i in 1:3){    #which temp variable
  for(j in 0:14){    #which offset
    Te=  ofset(Temps[[i]], j)[,0:85] 
    Pairs = rbind(Te[1,fn1(MF[1,]):ncol(Te)],MF[1,fn1(MF[1,]):ncol(Te)])
    for(k in 2:27){
      M = rbind(Te[k,fn1(MF[k,]):ncol(Te)],
                MF[k,fn1(MF[k,]):ncol(Te)])
      Pairs = cbind(Pairs,M)}
    Pairs = t(Pairs)
    #LM
    line = lm(Pairs[,2]~Pairs[,1])$coefficients
    Intercepts[[l]][i,j+1] = line[1]
    Slopes[[l]][i,j+1] = line[2]
    ell = nrow(Pairs)
    ##GRAPHING
    plot(Te[1,fn1(MF[1,]):ncol(Te)], MF[1,fn1(MF[1,]):ncol(Te)],
         xlim = c(-5,110), ylim = c(min(MF),max(MF)) ,
         main = paste('M.F. v.', Tnames[i]),
         xlab = paste(Tnames[i], 'offset by', j),
         ylab = "Multiplicative Factor",pch = 20)
    for(k in 2:27){
      points(Te[k,fn1(MF[k,]):ncol(Te)],MF[k,fn1(MF[k,]):ncol(Te)], col = colors[k], pch = 20)
    }
    lines(x,line[1] + line[2]*x, col = "red")
    ##TESTING
    test = cor.test(Pairs[,2], Pairs[,1],
                    method = whichmethod[l], exact = TRUE)
    p = round(test$p.value,5)
    rho = round(test$estimate,5)
    rhos[[l]][i,j+1] = rho
    pvals[[l]][i,j+1] = p
    print(paste(Tnames[i], 'offset by', j,': p-value = ', p, 
                'cor = ', rho, whichmethod[l], ell, "Pairs"))
  }
}}

par(mfrow = c(2,3), mar = c(4,4,4,4))
for(l in 1:3){
for(i in 1:3){
  plot(1:15,rhos[[l]][i,1:15], pch = 19,
       ylim = c(-0.2,0.2),
       main = c(paste( whichmethod[l],"correlation"),
       paste("for M.F.(t) v.",Tnames[i])),
       ylab = "Correlation", xlab = "Offset in Days")
  abline(h = 0, col = "grey")
}
for(i in 1:3){
  plot(1:15,pvals[[l]][i,1:15], pch = 19,
       main = paste("p-value v. offset"),
       ylab = "p-value", xlab = "Offset in Days")
  if(BH(pvals[[l]][i,1:15], 0.05)[1] > 0){
    points(BH(pvals[[l]][i,1:15], 0.05),
           pvals[[l]][i,BH(pvals[[l]][i,1:15], 0.05)],
           col = "red", pch = 19)
  }
  abline(h = 0.05, col = "red")
}}


##LINEAR REGRESSION BETA
for(i in 1:3){
  plot(1:15,Intercepts[[l]][i,1:15], pch = 19,
       ylim = c(min(Intercepts[[i]]),max(Intercepts[[i]])),
       main = paste("Regression intercept v. Offset"),
       ylab = "Intercept", xlab = paste(Tnames[i],"Offset in Days"))
  abline(h = 0, col = "grey")
}

for(i in 1:3){
  plot(1:15,Slopes[[l]][i,1:15], pch = 19,
       ylim = c(min(Slopes[[i]]),max(Slopes[[i]])),
       main = paste("Regression Slope v. Offset"),
       ylab = "Slope", xlab = paste(Tnames[i],"Offset in Days"))
  abline(h = 0, col = "grey")
}


###########################################################################
#############################################################################
###########################################################################
#########################################################################
#############################DERIVATIVES################################

time = c(1:length(Places[1,]))   #x axis
J = cbind(time,t(Places[1,]))    #data matrix
dim(deriv(J,0))
DV = matrix(0, nrow = 27,ncol = 86)   #all the derivative data
for(i in 1:27){
  J = cbind(time,t(Places[i,]))
  DV[i,] = deriv(J,0)[,2]
}
dim(DV)


par(mfrow = c(2,2),mar = c(4,4,4,4))  #use DV > 0 since some places had one guy in quarantine (LA, cook)
Z = matrix(0,nrow = 3, ncol = 15)     #ncol = 15 for usual
rhos = list(Z,Z,Z)
pvals = list(Z,Z,Z)
Intercepts = list(Z,Z,Z)
Slopes = list(Z,Z,Z)
for(l in 1:3){
for(i in 1:3){
  for(j in 0:14){       #0:14 or 0:20
    Te=  ofset(Temps[[i]], j)#[,1:72]   
    Pairs = rbind(Te[1,DV[1,]>0],DV[1,DV[1,]>0])
    for(k in 2:27){
      M = rbind(Te[k,DV[k,]>0],
                DV[k,DV[k,]>0])
      Pairs = cbind(Pairs,M)
    }
    Pairs = t(Pairs)
    #LM COEFFICIENTS
    line = lm(Pairs[,2]~Pairs[,1])$coefficients
    #cpairs = mcenter(Pairs)
    #linec = lm(cpairs[,2]~cpairs[,1])$coefficients
    Intercepts[[l]][i,j+1] = line[1]
    Slopes[[l]][i,j+1] = line[2]
    ell = nrow(Pairs)
    # #GRAPHING
    plot(Te[1,fnz(DV[1,]):ncol(Te)], DV[1,fnz(DV[1,]):ncol(Te)],
     xlim = c(0,110), ylim = c(min(DV),max(DV)) ,
     main = paste('Derivative v.', Tnames[i]),
     xlab = paste(Tnames[i], 'offset by', j),
     ylab = "Derivative",pch = 20)
    for(k in 2:27){
    points(Te[k,fnz(DV[k,]):ncol(Te)],DV[k,fnz(DV[k,]):ncol(Te)], 
           col = colors[k], pch = 20)
    }
    lines(x,line[1] + line[2]*x, col = "red")
    #TESTING
    test = cor.test(Pairs[,2], Pairs[,1],
                   method = whichmethod[l], exact = TRUE)
    p = round(test$p.value,5)
    rho = round(test$estimate,5)
    rhos[[l]][i,j+1] = rho
    pvals[[l]][i,j+1] = p
    print(paste(Tnames[i], 'offset by', j,': p-value = ', p, 
                'cor = ', rho, whichmethod[l], ell, "Pairs"
    ))
  }
}}
mcenter = function(matrix){
  M = c(mean(matrix[,1]),mean(matrix[,2]))
  Out = matrix - matrix
  for(i in 1:nrow(matrix)){
    Out[i,] = matrix[i,] - M
  }
  return(Out)
}


par(mfrow = c(2,3), mar = c(4,4,4,4))
for(l in 1:3){
  for(i in 1:3){
    plot(1:15,rhos[[l]][i,1:15], pch = 19,
         ylim = c(-0.2,0.2),
         main = c(paste( whichmethod[l],"correlation"),
                  paste("for N'(t) v.",Tnames[i])),
         ylab = "Correlation", xlab = "Offset in Days")
    abline(h = 0, col = "grey")
  }
  for(i in 1:3){
    plot(1:15,pvals[[l]][i,1:15], pch = 19,
         main = paste("p-value v. offset"),
         ylab = "p-value", xlab = "Offset in Days")
    if(BH(pvals[[l]][i,1:15], 0.05)[1] > 0){
      points(BH(pvals[[l]][i,1:15], 0.05),
             pvals[[l]][i,BH(pvals[[l]][i,1:15], 0.05)],
             col = "red", pch = 19)
    }
    abline(h = 0.05, col = "red")
  }}


##LINEAR REGRESSION BETA
for(i in 1:3){
  plot(1:15,Intercepts[[l]][i,1:15], pch = 19,
       ylim = c(min(Intercepts[[i]]),max(Intercepts[[i]])),
       main = paste("Regression intercept v. Offset"),
       ylab = "Intercept", xlab = paste(Tnames[i],"Offset in Days"))
  abline(h = 0, col = "grey")
}

for(i in 1:3){
  plot(1:15,Slopes[[l]][i,1:15], pch = 19,
       ylim = c(min(Slopes[[i]]),0.5),
       main = paste("Regression Slope v. Offset"),
       ylab = "Slope", xlab = paste(Tnames[i],"Offset in Days"))
  abline(h = 0, col = "grey")
}




###########################################################################
#############################################################################
###########################################################################
#########################################################################
#############################################################################
##########################################################################'
#Which places do  we want to look at? 


dim(Data)
par(mfrow=c(1,1))
stot = Data[,97]            #case totals
length(stot[stot > 3000])   #how many counties have >3000 cases

for(i in 1:3253){           #let's see which ones
  if(stot[i] > 2000){
    print(Data[i,6:7])
    print(stot[i])
  }
}

plot(1:3255,Data[,92])



#####################################################################
########LOCATION : LAT/LONG##########################################
#####################################################################
####################################################################
#####################################################################
####################################################################
#####################################################################
####################################################################
#
Data2 = Data[, 9:ncol(Data)]
Lat = Data2[,1]

Data2 = Data2[Lat > 25,]
Lat = Data2[,1]

Data2 = Data2[Lat < 50,]      #continental US
Lat = Data2[,1]
Long = Data2[,2]


par(mfrow = c(1,1), mar = c(4,4,4,4))
plot(Long, Lat, pch = 20)


TS = Data2[,-c(1,2,3)]
dim(TS)


time = c(1:length(TS[1,]))   #x axis
J = cbind(time,t(TS[1,]))    #data matrix
dim(deriv(J,0))
DV = matrix(0, nrow = 3110,ncol = 81)   #all the derivative data
for(i in 1:3110){
  J = cbind(time,t(TS[i,]))
  DV[i,] = deriv(J,0)[,2]
}


Y = TS[,81]
plot(Lat,Y) 

dim(TS)
dim(DV)

rhos = c()
pvals = c()
for(i in 1:81){
  test = cor.test(Lat, DV[,i],
           method = "kendall", exact = TRUE)
  p = round(test$p.value,10)
  rho = round(test$estimate,10)
  rhos[i] = rho
  pvals[i] = p
  print(paste('day ', i, 
              'rho = ', rho,'p-value = ', p
  ))
}


par(mfrow = c(2,1))
plot(time,rhos, pch = 19)
abline(h = 0, col = "red")
plot(time,pvals, pch = 19)
if(BH(pvals, 0.025)[1] > 0){
  points(BH(pvals, 0.025), pvals[BH(pvals, 0.025)],
         col = "red", pch = 20)
}
abline(h = 0.05, col = "red")






par(mfrow = c(2,3), mar = c(4,3,4,2))
for(i in 1:3){
  plot(1:15,rhos[i,1:15], pch = 19,
       main = paste("Rho for N'(t) &",Tnames[i]),
       ylab = "Rho", xlab = "Offset in Days")
  abline(h = 0, col = "red")
}
for(i in 1:3){
  plot(1:15,pvals[i,1:15], pch = 19,
       main = paste("p-value for N'(t) &",Tnames[i]),
       ylab = "p-value", xlab = "Offset in Days")
  if(BH(pvals[i,1:15], 0.05)[1] > 0){
    points(BH(pvals[i,1:15], 0.05), pvals[i,BH(pvals[i,1:15], 0.05)],
           col = "red", pch = 19)
  }
  abline(h = 0.05, col = "red")
}

