
% Generation if Figures similar to Fig3,4 in [1]

%INPUT: 
%load optimal fits
load('BestFittingTauMin0I2Std3April8.mat')
%Load NSB entropy estimates for the optimal completion distributions
load('NSBEntropyDataI2Std3.mat')
%load('BayesianSelImportSamplingFeb10.mat')

C=5;
nb=11;
MeanCurves=cell(1,C);
ErrorBarCurves=cell(1,C);
LogIntTot=zeros(1,C)

%Load curve estimates of the PDF and its respective errorbars, estimated
%using importance sampling
for i=1:C
    i
  Name4=['ImportSampCorr2AndBayesSelPurkinjeI2Std3April8C'  num2str(i) '.mat'];
  load(Name4)
  MeanCurves{1,i}=CurveMean;
  ErrorBarCurves{1,i}=sqrt(ErrorBarCurve);
  LogIntTot(i)=LogInt;
  
end


fig=figure
red_level = 0.5;
blue_level = 0.5;
green_level = 0.5;
cyan_level = 0.5;
% yellow_level = 0.5;
magenta_level = 0.5;
purple_level= 0.5;



shift{1} = [red_level 0 0];
shift{2} = [0 0 blue_level];
shift{3} = [0 green_level 0];
shift{5} = [0 cyan_level cyan_level];
shift{4} = [purple_level 0 purple_level];
shift{6} = [0.8 0.6 0.2];
shift{7} = [0.4 0.5 0.6];


shiftStair_color = [magenta_level 0 magenta_level];
baseline_color = [0.5 0.2 0.1];
 
capsize = 0;
handaxes1 = axes('position', [0.1 0.1 0.8 0.8]);
set(gcf,'color','w');

PD=nf/sum(nf);
EB=sqrt(nf)/sum(nf);



shadeUpEB = movmean(PD+EB,nb);

shadeDnEB= movmean(PD-EB,nb);


p(1)=plot(UTimes,PD,'.k');
%p(1)=errorbar(UTimes,PD,EB,'.k') ;
legendInfo=cell(1,5);
legendInfo{1} = ['Data, $H_{0}=$' num2str(S_nsb, ' %10.2f') '$\pm$ ' num2str(dS_nsb, ' %10.2f') ];
hold on
fill([UTimes,fliplr(UTimes)],[shadeUpEB,fliplr(shadeDnEB)], 'k','FaceAlpha',0.3,'EdgeColor','none');

v=[1,2,3,5,4];
%v=[1,2,4,5,3];

for l=1:5
  i=v(l)
%[ProbMod]= ModelMultiGammaDetermDt(XBest{1,j},j,UTimes,Dt);
hold on
p(i+1)=plot(UTimes,MeanCurves{1,i},'Color',shift{i},'LineWidth',2.5);
shadeUp = MeanCurves{1,i}+ErrorBarCurves{1,i};
shadeDn= MeanCurves{1,i}-ErrorBarCurves{1,i};
hold on
num2str(i)
A=sprintfc( ', $H_{%d} = $', i);
B=char(A)
legendInfo{i+1} = ['$M = $' num2str(i) B num2str(-LogLikelihood(i)/sum(nf), ' %10.2f')];
fill([UTimes,fliplr(UTimes)],[shadeUp,fliplr(shadeDn)],shift{i},'FaceAlpha',0.3,'EdgeColor','none');
clear shadeUp shadeDn
end

h=legend([p(1) p(2) p(3) p(4) p(5) p(6)], legendInfo)
set(h,'Interpreter','latex')
rect = [0.68, 0.27, .20, .2];
set(h, 'Position', rect)


ax=gca;
ax.LineWidth=2;
ax.XLim = [min(UTimes) max(UTimes)]; %Fix this according to Utimes ranges
ax.YAxis.Exponent = -3;
%set(gca,'XScale','log');
%set(gca,'YScale','log');
set(gca, 'box', 'off');
%title('$I=0.1$ nA','FontSize',52,'Interpreter','latex')
xlabel('Inter-spike interval (ms)','FontSize',38,'Interpreter','latex')
ylabel('PDF','FontSize',38,'Interpreter','latex')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 35)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 35);



handaxes2 = axes('position', [0.5 0.6 0.35 0.35]);
y=nf/sum(nf);
Dy=sqrt(nf)/sum(nf);

%ErrorNeg=log(y)-log(y-Dy);
%Indd=(ErrorNeg==Inf);
%ErrorNeg(Indd)=8; %Change this accordingly with every Concentration!!
%ErrorPos=log(y+Dy)-log(y);
%errorbar((UTimes),log(y),ErrorNeg,ErrorPos,'.k');
plot((UTimes),log10(y),'.k')
hold on
UpE=movmean(log10(y+Dy),nb);
DownE=movmean(log10(y-Dy),nb);
Indd=(DownE==-Inf);
DownE(Indd)= -18;
%DownE=movmean(DownE,23);

fill([UTimes,fliplr(UTimes)],[UpE,fliplr(DownE)],'k','FaceAlpha',0.3,'EdgeColor','none');

for l=1:5
 i=v(l);   
%[ProbMod]= ModelMultiGammaDetermDt(XBest{1,j},j,UTimes,Dt);
yt=MeanCurves{1,i};
Dyt=ErrorBarCurves{1,i};
hold on
plot((UTimes),log10(yt),'Color',shift{i},'LineWidth',2.5);
shadeUp = log10(yt+Dyt);
shadeDn= log10(yt-Dyt);
fill([UTimes,fliplr(UTimes)],[shadeUp,fliplr(shadeDn)],shift{i},'FaceAlpha',0.3,'EdgeColor','none');
clear shadeUp shadeDn yt Dyt
end
ax = gca;
ax.LineWidth=1.5;
ax.YLim = [-7 -1];
ax.XLim = [min(UTimes) max(UTimes)];
set(gca,'XScale','log');
%set(gca,'YScale','log');
xlabel('Inter-spike interval (ms)','FontSize',16,'Interpreter','latex')
ylabel('Log(PDF)','FontSize',26,'Interpreter','latex')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 26)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 26)
yticks([-7 -5 -3 -1])
yticklabels({'-7','-5','-3','-1'})
orient(fig,'landscape')
%print(fig,'PaperFigs/BestFitsPurkinjeI0p1Std3Feb2019nb11.eps','-depsc2')
