

% Use this file for those models that have a trivial path (Ex. C=5 for Purkinje  )
Name2=['BestFittingForFigNTauMin0I0p7Std3C5Sept2019N92915.mat'];
load(Name2)

k=4;
Name=['HessianI0p7Std3C' num2str(k) 'MathemForFigNSept2019N92915.dat']
Hess=load(Name);
CovN=inv(-Hess);
CovN=(CovN+CovN.')/2;

%Generate random number:
Fmax=-LogObjectFuncMultiGammaModelDeterm(XBest{1,k},UTimes,k,nf,Dt,TauMin,Z1,Z2);
N=1e6;
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
 else
 
 j=j+1;
 %These expressions had an error
 end
 
end
LogInt= log(int)+Fmax-log(N);

LogStdError=sqrt( VarSum/N -(int/N)^2)/(int/sqrt(N));

CurveMean=CurveMeanT/int;
ErrorBarCurve=(CurveErrT1/int) -CurveMean.^2; 
%ErrorBarCurve= (CurveErrT1 - 2*CurveMean.*CurveErrT2+ (CurveMean.^2).*CurveErrT3)/(int^2);
filename=['ImportSamplesPurkinjeI0p7Std3C' num2str(k) 'Sept2019N92915.mat'];
save(filename,'RSamp')
filename2=['ImportSampAndBayesSelPurkinjeI0p7Std3C' num2str(k) 'Sept2019N92915.mat'];
save(filename2,'LogInt','LogStdError','CurveMean','ErrorBarCurve','s','j','N')


