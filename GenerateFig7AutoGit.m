% Generate figure 7 in paper. This file is a copy of FigureHistJSD.m

% Two separate data sets. Bootstrapping one of them.
% Number of bootstraps
NB= 100;
% Loading data set 2
load('DataPurkinjeI2p5Std3April8V2.mat')

UTimes2=UTimes;
nf2=nf;
SpikeTimes2=SpikeTimes;
clear UTimes nf SilencesN SpikeTimes

% Loading data set 1
load('BestFittingTauMin0I2p5Std3April8.mat')

PData=nf/sum(nf);
[ProbFit]= ModelMultiGammaDetermDt(XBest{1,4},4,UTimes,Dt);

% Boostrapping with replacement from Data set 2 

for i=1:NB
Boots = datasample(SpikeTimes2,length(SpikeTimes2));
MIN=min(min(UTimes),min(Boots));
MAX=max(max(UTimes),max(Boots));
edges=[MIN:Dt:MAX];
hT1 = hist(SpikeTimes,edges);

hT2 = hist(Boots,edges);
PData1=hT1/sum(hT1);
PData2=hT2/sum(hT2);
PM=(PData1+PData2)/2;

% JSD (Data 1 | Boots Data 2)
JSData(i)=sum(log(PData1.^(PData1))-log(PM.^(PData1)))/2 +sum(log(PData2.^(PData2))-log(PM.^(PData2)))/2
JSDataElem=(log(PData1.^(PData1))-log(PM.^(PData1)))/2 +(log(PData2.^(PData2))-log(PM.^(PData2)))/2

%[ProbFit]= ModelMultiGammaDetermDt(XBest{1,4},4,edges,Dt);
[XFitNew ChiFitBoot, I, xT, ChiT]=FittingHierarMultiGammaModelDeterm(edges,4,XBest{1,3},XBest,hT2,Dt,TauMin,Z1,Z2);
[ProbFitNew]= ModelMultiGammaDetermDt(XFitNew,4,edges,Dt);
PMNew=(PData2+ProbFitNew)/2;
% JSD (Fit Boots | Boots Data 2)
JSRealFitNew(i)=sum(log(PData2.^(PData2))-log(PMNew.^(PData2)))/2+ sum(log(ProbFitNew.^(ProbFitNew))-log(PMNew.^(ProbFitNew)))/2

load('XInterCurrent2p5Nov2019.mat')
[ProbPred]= ModelMultiGammaDetermDt(Xinter,4,edges,Dt);
PM3=(PData2+ProbPred)/2;
% JSD (Pred | Boots Data 2)
JSRealPred(i)=sum(log(PData2.^(PData2))-log(PM3.^(PData2)))/2+ sum(log(ProbPred.^(ProbPred))-log(PM3.^(ProbPred)))/2
%clear Boots hT2 PData2 PMNew PM3
end

% Generate plot of the JSD histograms

r=JSData,
r2=JSRealFitNew 
r3=JSRealPred

close all;
clc;

[f,xi] = ksdensity(r);
ht = histogram(r,15);
z = sum(ht.Values*ht.BinWidth);
stem(ht.BinEdges(1:end-1)+ht.BinWidth/2,ht.Values/z);
hold on;
p(1)=plot(xi,f,'Color','k','LineWidth',2.5);


hold on

[f2,x2] = ksdensity(r2);
ht2 = histogram(r2,15);
z2 = sum(ht2.Values*ht2.BinWidth);
stem(ht2.BinEdges(1:end-1)+ht2.BinWidth/2,ht2.Values/z2);
hold on
p(2)=plot(x2,f2,'Color','b','LineWidth',2.5);

hold on

[f3,x3] = ksdensity(r3);
ht3 = histogram(r3,15);
z3 = sum(ht2.Values*ht2.BinWidth);
stem(ht3.BinEdges(1:end-1)+ht3.BinWidth/2,ht3.Values/z2);
hold on
p(3)=plot(x3,f3,'Color','r','LineWidth',2.5);
%h=legend([p(1) p(2) p(3) p(4) p(5) p(6)], legendInfo)
ax=gca;
ax.LineWidth=2;
legendInfo{1}=['JSD(D $\|$ B)'];
legendInfo{2}=['JSD(F $\|$ B)'];
legendInfo{3}=['JSD(P $\|$ B)'];



h=legend([p(1) p(2) p(3)],legendInfo)
set(gcf,'color','w');
ax=gca;
ax.LineWidth=2;
ax.YAxis.Exponent = 2;
ax.XLim = [0.008 0.016];
ax.XTick = [0.008, 0.01,0.012, 0.014,0.016,0.018];
ax.XTickLabel = {'0.008','0.01', '0.012', '0.014','0.016','0.018'};
%ax.YTick = [0, 500,1000, 1500 ];
%ax.YTickLabel = {'0','500', '1000', '1500'};


xlabel('JSD','FontSize',48,'Interpreter','latex')
set(gca, 'FontSize', 48)
ylabel('PDF','FontSize',48,'Interpreter','latex')
set(gca, 'FontSize', 48)
title('$I=2.5$ nA','FontSize',58,'Interpreter','latex')
set(h,'Interpreter','latex')