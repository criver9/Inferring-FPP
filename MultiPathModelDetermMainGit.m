
%Finds optimal fits for the first M models in the multi-path model family

%INPUT:
%Data set
load('DataPurkinjeI3p3Std3April8.mat')
% Time resolution of data
Dt=0.025;
% Maximum model in the hierarchy to fit
M=3;

% Parameters of the Prior distributions
TauMin=0;
ZTau=20;
ZL=20;


% Propose starting points for optimization. Can be previous solutions of
% this functions. 
load('XPrevI3p3April8.mat');

%OUTPUT

% Model parameter for optimal fits for each of the first M models: XBest
% Log-Likelihood for the fitted models: LogLikelihood


%------------------------------------
XBest=cell(1,M);
ChiBest=zeros(1,M);
LogLikelihood=zeros(1,M);

Xprev=[];


for c=1:M
   
   [XBest{1,c} ChiBest(c), I, xT{1,c}, ChiT{1,c}]=FittingHierarMultiGammaModelDeterm(UTimes,c,Xprev,XPrev2,nf,Dt,TauMin,ZTau,ZL)
   [LogLikelihood(c)]= LogLikelihoodNewLoss(XBest{1,c},UTimes,c,nf,Dt)
   Xprev=XBest{1,c};
   
end


%save('BestFittingTauMin0I1Std3C5April8')
