%For the Best Generalized Bayesian model selection, plots parameters as a
%function of Current. BEst model has 4 paths.

Best=4; % From BayesianSelGenWMathemFeb7.mat
In1=1:3:3*Best-1;
In2=2:3:3*Best-1;
In3=3:3:3*Best-1;
Current=[0.1, 0.5, 0.7, 1, 2,3];
BestParam=zeros(length(Current),3*Best-1);
Name{1}='BestFittingTauMin0I0p1Std3Feb11.mat';
Name{2}='BestFittingTauMin0I0p5Std3Feb11.mat';
Name{3}='BestFittingTauMin0I0p7Std3April8.mat';
Name{4}='BestFittingTauMin0I1Std3April8.mat';
Name{5}='BestFittingTauMin0I2Std3April8.mat';
Name{6}='BestFittingTauMin0I3Std3April8.mat';

Name2{1}='HessianI0p1Std3C4MathemFeb11.dat'
Name2{2}='HessianI0p5Std3C4MathemFeb11.dat'
Name2{3}='HessianI0p7Std3C4MathemApril8.dat'
Name2{4}='HessianI1Std3C4MathemApril8.dat'
Name2{5}='HessianI2Std3C4MathemApril8.dat'
Name2{6}='HessianI3Std3C4MathemApril8.dat'


for l=1:6
    l
    
    load(Name{l})
    
    Hess=load(Name2{l});
    XT=XBest{1,Best};
    BestParam(l,:)=XT;
    TAveTemp=BestParam(l,In1).*BestParam(l,In2);
    [M,I]=sort( TAveTemp);
    I
    ParamTauTemp=BestParam(l,In1);
    ParamTau(l,:)=ParamTauTemp(I);
    ParamLTemp=BestParam(l,In2);
    ParamL(l,:)=ParamLTemp(I);
    TAve(l,:)=TAveTemp(I); 
    CoefVarSq(l,:)=1./ParamL(l,:);
    VarTAveTemp=(BestParam(l,In1).^2).*BestParam(l,In2);
    VarTAve(l,:)=VarTAveTemp(I);
    ParamFluxTemp=[1 BestParam(l,In3)];
    ParamFlux(l,:)=ParamFluxTemp(I);
    ParamProb(l,:)=ParamFlux(l,:)./(sum(ParamFluxTemp));
end
 
%% . Extrapolation:
Conc=[0.1 0.5 0.7 1 2 3];
ti=3.5;
TAve
TAveInter(1) = interp1(Conc',TAve(:,1),ti,'linear','extrap');
TAveInter(2) = interp1(Conc',TAve(:,2),ti,'linear','extrap');
TAveInter(3) = interp1(Conc',TAve(:,3),ti,'linear','extrap');
TAveInter(4) = interp1(Conc',TAve(:,4),ti,'linear','extrap');


CVInter(1) = interp1(Conc',CoefVarSq(:,1),ti,'linear','extrap');
CVInter(2) = interp1(Conc',CoefVarSq(:,2),ti,'linear','extrap');
CVInter(3) = interp1(Conc',CoefVarSq(:,3),ti,'linear','extrap');
CVInter(4) = interp1(Conc',CoefVarSq(:,4),ti,'linear','extrap');

ProbInter(1) = interp1(Conc',ParamProb(:,1),ti,'linear','extrap');
ProbInter(2) = interp1(Conc',ParamProb(:,2),ti,'linear','extrap');
ProbInter(3) = interp1(Conc',ParamProb(:,3),ti,'linear','extrap');
ProbInter(4) = interp1(Conc',ParamProb(:,4),ti,'linear','extrap');

%% Extrapolation parameter vector

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
%% PLot the PDF data, Fitted and Predicted
load('DataPurkinjeI3p5Std3April8.mat')
Dt=0.025;
[ProbCurve]= ModelMultiGammaDetermDt(Xinter,4,UTimes,Dt)
plot(UTimes,nf/sum(nf),'.')
hold on
plot(UTimes,ProbCurve,'r')
load('BestFittingTauMin0I3p5Std3April8.mat')
[Probfit]= ModelMultiGammaDetermDt(XBest{1,4},4,UTimes,Dt)
plot(UTimes,Probfit,'k')
legend('Data','Predicted','Fitted')
%% Kl Divergence
TT=Dt:Dt:1000;
[ProbPr]= ModelMultiGammaDetermDt(Xinter,4,UTimes,Dt);
[Probfit]= ModelMultiGammaDetermDt(XBest{1,4},4,UTimes,Dt);
KL=sum( Probfit.*log(Probfit./ProbPr));



%% Plot Paths together 
Colors{1} = [0 0 0.5];
Colors{2} = [0.8 0.6 0.2];
Colors{3} = [0 0.5 0.5];
Colors{4} = [0.5 0 0.2];
Current=[0.1, 0.5, 0.7, 1, 2,3, 2.5];

fig=figure;

pos1 = [0.07 0.4 0.4 0.4];
subplot(1,3,1)

h=plot(Current,[TAve(:,1);TAveInter(1)]  ,'Color',Colors{1},'LineWidth',1);
hold on
plot(Current,[TAve(:,2);TAveInter(2)],'Color',Colors{2},'LineWidth',1);
hold on
plot(Current,[TAve(:,3);TAveInter(3)],'Color',Colors{3},'LineWidth',1);
hold on
plot(Current,[TAve(:,4);TAveInter(4)],'Color',Colors{4},'LineWidth',1);

title({'Completion Time '},'FontSize',14,'FontName','Arial')
xL=xlabel({'Current (nA)'},'FontSize',12,'FontName','Arial')
yL=ylabel({'$\overline{T}=\tau L$'},'Interpreter','latex')

ax = ancestor(h, 'axes');
xrule = ax.XAxis;
ax.XTick = [0.1, 0.5,0.7, 1, 2,3];
ax.XTickLabel = {'0.1','', '0.7', '1', '2','3'};
xrule.FontSize = 10;
xL.FontSize = 12;
yrule=ax.YAxis
yrule.FontSize=12;
yL.FontSize=14;


set(gcf,'color','w');
axis square

pos2 = [0.4 0.4 0.3 0.3];
subplot(1,3,2)
%subplot('Position',pos2)
h=plot(Current,[CoefVarSq(:,1);CVInter(1) ],'Color',Colors{1},'LineWidth',1);
hold on
plot(Current,[CoefVarSq(:,2);CVInter(2)],'Color',Colors{2},'LineWidth',1);
hold on
plot(Current,[CoefVarSq(:,3);CVInter(3)],'Color',Colors{3},'LineWidth',1);
hold on
plot(Current,[CoefVarSq(:,4);CVInter(4)],'Color',Colors{4},'LineWidth',1);
%title('VarTAve')
title({'Coefficient of Variation^2 '},'FontSize',14,'FontName','Arial')
xL=xlabel('Current (nA)','FontSize',12,'FontName','Arial')
yL=ylabel({'$CV^2=1/L$'},'Interpreter','latex')

ax = ancestor(h, 'axes');

xrule = ax.XAxis;
ax.XTick = [0.1, 0.5,0.7, 1, 2,3];
ax.XTickLabel = {'0.1','', '0.7', '1', '2','3'};
xrule.FontSize = 10;
xL.FontSize = 12;
yrule=ax.YAxis
yrule.FontSize=12;
yL.FontSize=14;


set(gcf,'color','w');
axis square
Current=[0.1, 0.5, 0.7, 1, 2,3];
pos3 = [0.7 0.4 0.3 0.3];
%subplot('Position',pos3)
subplot(1,3,3)
h=plot(Current,ParamProb(:,1),'Color',Colors{1},'LineWidth',1);
hold on
plot(Current,ParamProb(:,2),'Color',Colors{2},'LineWidth',1);
hold on
plot(Current,ParamProb(:,3),'Color',Colors{3},'LineWidth',1);
hold on
plot(Current,ParamProb(:,4),'Color',Colors{4},'LineWidth',1);
%title('Prob')
title({'Probability of path'},'FontSize',14,'FontName','Arial')
xL=xlabel('Current (nA)','FontSize',12,'FontName','Arial')
yL=ylabel('$p$','Interpreter','latex','FontSize',16,'FontName','Arial')
ax = ancestor(h, 'axes');

xrule = ax.XAxis;
ax.XTick = [0.1, 0.5,0.7, 1, 2,3];
ax.XTickLabel = {'0.1','', '0.7', '1', '2','3'};
xrule.FontSize = 10;
xL.FontSize = 14;
yrule=ax.YAxis
yrule.FontSize=12;
yL.FontSize=14;
axis square
