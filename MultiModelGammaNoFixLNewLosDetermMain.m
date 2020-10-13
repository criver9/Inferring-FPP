%load('../../NEURON SIMULATOR/prknj/DataPurkinjeI0p7Std3Sept2019N1e6.mat')
load('DataPurkinjeI3p3Std3April8.mat')
Dt=0.025;
C=5;
XBest=cell(1,C);
ChiBest=zeros(1,C);
LogLikelihood=zeros(1,C);
%xprev2 is the best fitting from first lost function.
Xprev=[];
load('XPrevI3p3April8.mat');
%Xprev=XPrev2{1,4};
TauMin=0;
Z1=20;
Z2=20;


for c=1:5
    c
    
   [XBest{1,c} ChiBest(c), I, xT{1,c}, ChiT{1,c}]=FittingHierarMultiGammaModelDeterm(UTimes,c,Xprev,XPrev2,nf,Dt,TauMin,Z1,Z2)
   %[XBest{1,c} ChiBest(c), I, xT{1,c}, ChiT{1,c}]=FittingHierarMultiGammaModelDetermfmincon(UTimes,c,Xprev,XPrev2,nf,Dt,TauMin,Z1,Z2)
   [LogLikelihood(c)]= LogLikelihoodNewLoss(XBest{1,c},UTimes,c,nf,Dt)
   Xprev=XBest{1,c};
    %save('TempFileI1Determ')
    
end

%save('BestFittingTauMin0I1Std3C5April8')

