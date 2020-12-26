%Estimates the logarithms of the marginal likelihood and the best fitted curves for the probability distribution for a given model in
%the multi-path model family. 

% Download the best fits for each model obtained using 
% MultiModelGammaNoFixLNewLosDetermMain.m. 
Name2=['BestFittingTauMin0I0p7Std3April8.mat'];
load(Name2)
% Specify the M for which you want to estimate the marginal likelihood
k=5;
Name=['HessianI0p7Std3C' num2str(k) 'MathemApril8.dat']
Hess=load(Name);
%Estimate Covariance matrix 
CovN=inv(-Hess);
CovN=(CovN+CovN.')/2;

% Fmax corresponds to the value of log(P(D|Theta*,M))+Log(P(Theta*|M)),
% where theta* are the optimal parameters estimated in
% MultiModelGammaNoFixLNewLosDetermMain.m.
Fmax=-LogObjectFuncMultiGammaModelDeterm(XBest{1,k},UTimes,k,nf,Dt,TauMin,Z1,Z2);
%Number of samples to estimate importance sampling
N=1e5;

% Read notes on Importance sampling to understand this notation
int=0;
VarSum=0;
CurveMeanT=zeros(1,length(UTimes));
CurveErrT1=zeros(1,length(UTimes));

s=0;
j=0;
RSamp=zeros(N,length(XBest{1,k}));
while s<=N
 R = mvnrnd(XBest{1,k},CovN,1); %Save this only for the best model selected
 if sum(R<0)==0
 s=s+1;    
 RSamp(s,:)=R;   
 PdfR= mvnpdf(R,XBest{1,k},CovN);
 F=-LogObjectFuncMultiGammaModelDeterm(R,UTimes,k,nf,Dt,TauMin,Z1,Z2);
 %exp(F-Fmax)
 %F=- (R*CovN*R')+Fmax;
 % int and VarSum are for Bayesian model selection
 int=int+exp(F-Fmax)/(PdfR);
 VarSum=VarSum+ exp(2*(F-Fmax))/(PdfR^2); %See notes Important Sampling March 29/2018
 % Mean and Error of Curves Posterior Sampling
 ProbMod = ModelMultiGammaDetermDt(R,k,UTimes,Dt);
 CurveMeanT = CurveMeanT + (ProbMod*exp(F-Fmax)/PdfR);
 CurveErrT1=CurveErrT1 + (ProbMod.^2*exp(F-Fmax)/PdfR);
 end
 j=j+1;
end
%Logarithm of the marginal likelihood
LogInt= log(int)+Fmax-log(N)+log(factorial(k-1))-log(j/N);
% Read notes to understand terminology
VarT = VarSum -int^2/N;
SigmaT=sqrt(VarT)
% upper bound of error bar of the logarithm of the marginal likelihood
DP=Fmax-log(N)+log(factorial(k-1))-log(j/N)+log(int+SigmaT);
% Lower bound of error bar of the logarithm of the marginal likelihood
DN=Fmax-log(N)+log(factorial(k-1))-log(j/N)+log(int-SigmaT);
% Mean of the fitted curved
CurveMean=CurveMeanT/int;
%Error bar on the mean fitted curve
ErrorBarCurve=(CurveErrT1/int) -CurveMean.^2; 
%save(filename,'RSamp')
filename2=['ImportSampCorr2AndBayesSelPurkinjeI0p7Std3April8C' num2str(k) '.Dic2020N1e5.mat'];
save(filename2,'LogInt','DP','DN','CurveMean','ErrorBarCurve','j','s')