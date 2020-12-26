% GEneration of predicted probability distributions by interpolation of
% parameters

% INTERPOLATION:

% load samples from posterior distribution
Samples=cell(1,6);
Name=['/Users/catalina/Library/Mobile Documents/com~apple~CloudDocs/Purkinje Posterior Samples/'];
Name1=[Name 'ImportSamplesCorrPurkinjeI0p1Std3Feb11C4']
load(Name1)
Samples{1,1}=RSamp;
clear RSamp;
Name2=[Name 'ImportSamplesCorrPurkinjeI0p5Std3Feb11C4']
load(Name2)
Samples{1,2}=RSamp;
clear RSamp;
Name3=[Name 'ImportSamplesCorrPurkinjeI0p7Std3April8C4']
load(Name3)
Samples{1,3}=RSamp;
clear RSamp;
Name4=[Name 'ImportSamplesCorrPurkinjeI1Std3April8C4']
load(Name4)
Samples{1,4}=RSamp;
clear RSamp;
Name5=[Name 'ImportSamplesCorrPurkinjeI2Std3April8C4']
load(Name5)
Samples{1,5}=RSamp;
clear RSamp;
Name6=[Name 'ImportSamplesCorrPurkinjeI3Std3April8C4']
load(Name6)
Samples{1,6}=RSamp;
clear RSamp;


In1=2:3:3*Best;
In2=3:3:3*Best;
In3=1:3:3*Best;

for K=1:6
RSamples=Samples{1,K};  
TAveN=RSamples(:,In1).*RSamples(:,In2);  % Completion Time per Sample
TAve(K,:)=mean(TAveN,1);                 % Mean of completion time
StdTAve(K,:)=std(TAveN,0,1);             % Std of Completion Time
CVN=1./RSamples(:,In2);                  % Coeff of Variation per Sample
CVMean=mean(CVN,1);                      % Mean of Coeff of Variation
CVStd=std(CVN,0,1);                      % Std of Coeff of Variation 
Prob=RSamples(:,In3);                       % Probability per path
MeanProb(K,:)=mean(Prob,1);               % Mean Prob
StdProb(K,:)=std(Prob,0,1);               % Std Prob
TAveQuant(:,:,K)=quantile(TAveN,[0.2 0.5 0.8],1); % Quantiles for Tave
CVQuant(:,:,K)=quantile(CVN,[0.2 0.5 0.8],1);     % Quantile for Coeff Var
ProbQuant(:,:,K)=quantile(Prob,[0.2 0.5 0.8],1);  % Quantile for Probability
%Y = quantile(Prob,20)
end
% . Extrapolation:
Conc=[0.1 0.5 0.7 1 2 3];
Ti=2.5;

TAve
TAveInter(1) = interp1(Conc',squeeze(TAveQuant(2,1,:)),Ti,'linear','extrap');
TAveInter(1) = interp1(Conc',squeeze(TAveQuant(2,1,:)),Ti,'linear');
TAveInter(2) = interp1(Conc',squeeze(TAveQuant(2,2,:)),Ti,'linear');
TAveInter(3) = interp1(Conc',squeeze(TAveQuant(2,3,:)),Ti,'linear');
TAveInter(4) = interp1(Conc',squeeze(TAveQuant(2,4,:)),Ti,'linear');


CVInter(1) = interp1(Conc',squeeze(CVQuant(2,1,:)),Ti,'linear');
CVInter(2) = interp1(Conc',squeeze(CVQuant(2,2,:)),Ti,'linear');
CVInter(3) = interp1(Conc',squeeze(CVQuant(2,3,:)),Ti,'linear');
CVInter(4) = interp1(Conc',squeeze(CVQuant(2,4,:)),Ti,'linear');

ProbInter(1) = interp1(Conc',squeeze(ProbQuant(2,1,:)),Ti,'linear');
ProbInter(2) = interp1(Conc',squeeze(ProbQuant(2,2,:)),Ti,'linear');
ProbInter(3) = interp1(Conc',squeeze(ProbQuant(2,3,:)),Ti,'linear');
ProbInter(4) = interp1(Conc',squeeze(ProbQuant(2,4,:)),Ti,'linear');

% Extrapolation parameter vector

LInter=1./CVInter;
TauInter=TAveInter./LInter;
ProbInter=ProbInter/sum(ProbInter);
ProbFluxInter=ProbInter./ProbInter(1);


Xinter(1)=TauInter(1);
Xinter(2)=LInter(1);
Xinter(3)=ProbFluxInter(2);
Xinter(4)=TauInter(2);
Xinter(5)=LInter(2);
Xinter(6)=ProbFluxInter(3);
Xinter(7)=TauInter(3);
Xinter(8)=LInter(3);
Xinter(9)=ProbFluxInter(4);
Xinter(10)=TauInter(4);
Xinter(11)=LInter(4);

%GENERATE PLOT FIG 6 OF PAPER

% Plot PDF of Data, fitted and predicted
load('DataPurkinjeI2p5Std3April8.mat')
Dt=0.025;
nb=11;
[ProbCurve]= ModelMultiGammaDetermDt(Xinter,4,UTimes,Dt)
PD=nf/sum(nf);
EB=sqrt(nf)/sum(nf);
shadeUpEB = movmean(PD+EB,nb);
shadeDnEB= movmean(PD-EB,nb);
p(1)=plot(UTimes,PD,'.k');
hold on
fill([UTimes,fliplr(UTimes)],[shadeUpEB,fliplr(shadeDnEB)], 'k','FaceAlpha',0.3,'EdgeColor','none');
hold on
p(2)=plot(UTimes,ProbCurve,'r','LineWidth',2)
load('BestFittingTauMin0I2p5Std3April8.mat')
[Probfit]= ModelMultiGammaDetermDt(XBest{1,4},4,UTimes,Dt)
hold on
p(3)=plot(UTimes,Probfit,'b','LineWidth',2)

ax=gca;
ax.LineWidth=2;
ax.XLim = [2 30]; %Fix this according to Utimes ranges
ax.YAxis.Exponent = -3;
%set(gca,'XScale','log');
set(gca, 'box', 'off');
xlabel('Inter-spike interval (ms)','FontSize',48)
ylabel('PDF','FontSize',48,'Interpreter','latex')
title('$I=2.5$ nA','FontSize',58,'Interpreter','latex')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 48)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 48);
set(gcf,'color','w');
legend('Data','Predicted','Fitted')
legendInfo{1} = ['Data'];
legendInfo{2} = ['Predicted']
legendInfo{3} = ['Fitted']
h=legend([p(1) p(2) p(3)], legendInfo)

