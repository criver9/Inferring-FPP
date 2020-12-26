

%Figure of Parameters vs Current, using Importance sampling to sample from
%the posterior to estimate parameters. 
%NOTE: Samples are already organized by TAve.


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

Best=4;
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
 

% Plot Paths together Using Quantile as errorbars

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

Current=[0.1, 0.5, 0.7, 1, 2,3];

fig=figure(1);
set(gcf,'color','w')
%pos1 = [0.07 0.4 0.4 0.4];
%subplot(1,3,1)
%subplot('Position',pos1)
h=errorbar(Current,squeeze(TAveQuant(2,1,:)),squeeze(TAveQuant(2,1,:)-TAveQuant(1,1,:)), squeeze(TAveQuant(3,1,:)-TAveQuant(2,1,:)),'Color',shift{1},'LineWidth',2);
hold on
errorbar(Current,squeeze(TAveQuant(2,2,:)),squeeze(TAveQuant(2,2,:)-TAveQuant(1,2,:)), squeeze(TAveQuant(3,2,:)-TAveQuant(2,2,:)),'Color',shift{2},'LineWidth',2);
hold on
errorbar(Current,squeeze(TAveQuant(2,3,:)),squeeze(TAveQuant(2,3,:)-TAveQuant(1,3,:)), squeeze(TAveQuant(3,3,:)-TAveQuant(2,3,:)),'Color',shift{3},'LineWidth',2);
hold on
errorbar(Current,squeeze(TAveQuant(2,4,:)),squeeze(TAveQuant(2,4,:)-TAveQuant(1,4,:)), squeeze(TAveQuant(3,4,:)-TAveQuant(2,4,:)),'Color',shift{6},'LineWidth',2);

%title('TAve')
title({'Completion Time '},'FontSize',52)
xL=xlabel('$I$ (nA)','FontSize',48,'Interpreter','latex')
yL=ylabel({'$\overline{T}=\tau L$'},'Interpreter','latex','FontSize',48)

ax=gca;
ax.LineWidth=2;
ax = ancestor(h, 'axes');
xrule = ax.XAxis;
ax.XLim = [0 3.05];
ax.XTick = [0.1, 0.5,0.7, 1, 2,3];
ax.XTickLabel = {'0.1','', '0.7', '1', '2','3'};
ax.YLim = [0 110];
ax.YTick = [0, 20,40,60,80,100];
ax.YTickLabel = {'0','20','40','60','80','100'};
xrule.FontSize = 48;
xL.FontSize = 48;
yrule=ax.YAxis
yrule.FontSize=48;
yL.FontSize=48;


set(gcf,'color','w');
%axis square

%pos2 = [0.4 0.4 0.3 0.3];
%subplot(1,3,2)
%subplot('Position',pos2)
figure(2)
h=errorbar(Current,squeeze(CVQuant(2,1,:)),squeeze(CVQuant(2,1,:)-CVQuant(1,1,:)),squeeze(CVQuant(3,1,:)-CVQuant(2,1,:)),'Color',shift{1},'LineWidth',2);
hold on
errorbar(Current,squeeze(CVQuant(2,2,:)),squeeze(CVQuant(2,2,:)-CVQuant(1,2,:)), squeeze(CVQuant(3,2,:)-CVQuant(2,2,:)),'Color',shift{2},'LineWidth',2);
hold on
errorbar(Current,squeeze(CVQuant(2,3,:)),squeeze(CVQuant(2,3,:)-CVQuant(1,3,:)), squeeze(CVQuant(3,3,:)-CVQuant(2,3,:)),'Color',shift{3},'LineWidth',2);
hold on
errorbar(Current,squeeze(CVQuant(2,4,:)),squeeze(CVQuant(2,4,:)-CVQuant(1,4,:)), squeeze(CVQuant(3,4,:)-CVQuant(2,4,:)),'Color',shift{6},'LineWidth',2);
%title('VarTAve')
title({'Squared Coefficient of Variation  '},'FontSize',52)
xL=xlabel('$I$ (nA)','FontSize',48,'Interpreter','latex')
yL=ylabel({'$CV^2=1/L$'},'Interpreter','latex')

ax=gca;
ax.LineWidth=2;
ax = ancestor(h, 'axes');
ax.XLim = [0 3.05];
xrule = ax.XAxis;
ax.XTick = [0.1, 0.5,0.7, 1, 2,3];
ax.XTickLabel = {'0.1','', '0.7', '1', '2','3'};
ax.YTick = [0, 0.2,0.4, 0.6];
ax.YTickLabel = {'0','0.2', '0.4', '0.6'};
ax.YLim = [0 0.65];
xrule.FontSize = 48;
xL.FontSize = 48;
yrule=ax.YAxis
yrule.FontSize=48;
yL.FontSize=48;


set(gcf,'color','w');
%axis square

%pos3 = [0.7 0.4 0.3 0.3];
%subplot('Position',pos3)
%subplot(1,3,3)
figure(3)

h=errorbar(Current,squeeze(ProbQuant(2,1,:)),squeeze(ProbQuant(2,1,:))-squeeze(ProbQuant(1,1,:)),squeeze(ProbQuant(3,1,:))-squeeze(ProbQuant(2,1,:)),'Color',shift{1},'LineWidth',2);
hold on
errorbar(Current,squeeze(ProbQuant(2,2,:)),squeeze(ProbQuant(2,2,:))-squeeze(ProbQuant(1,2,:)),squeeze(ProbQuant(3,2,:))-squeeze(ProbQuant(2,2,:)),'Color',shift{2},'LineWidth',2);
hold on
errorbar(Current,squeeze(ProbQuant(2,3,:)),squeeze(ProbQuant(2,3,:))-squeeze(ProbQuant(1,3,:)),squeeze(ProbQuant(3,3,:))-squeeze(ProbQuant(2,3,:)),'Color',shift{3},'LineWidth',2);
hold on
errorbar(Current,squeeze(ProbQuant(2,4,:)),squeeze(ProbQuant(2,4,:))-squeeze(ProbQuant(1,4,:)),squeeze(ProbQuant(3,4,:))-squeeze(ProbQuant(2,4,:)),'Color',shift{6},'LineWidth',2);
%title('Prob')
title({'Probability of path'},'FontSize',52)
xL=xlabel('$I$ (nA)','FontSize',48,'Interpreter','latex')
yL=ylabel('$p$','Interpreter','latex','FontSize',48)

ax=gca;
ax.LineWidth=2;
ax = ancestor(h, 'axes');
xrule = ax.XAxis;
ax.XLim = [0 3.05];
ax.XTick = [0.1, 0.5,0.7, 1, 2,3];
ax.XTickLabel = {'0.1','', '0.7', '1', '2','3'};
ax.YTick = [0, 0.2, 0.4, 0.6, 0.8];
ax.YTickLabel = {'0','0.2', '0.4', '0.6', '0.8'};
xrule.FontSize = 48;
xL.FontSize = 48;
yrule=ax.YAxis
yrule.FontSize=48;
yL.FontSize=48;
%axis square
set(gcf,'color','w');

%orient(fig,'landscape')
%print('PaperFigs/BestModelVsCurrent21Jan2019','-dpdf','-fillpage')












