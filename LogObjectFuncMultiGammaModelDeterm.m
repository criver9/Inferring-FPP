function [L]= LogObjectFuncMultiGammaModelDeterm(x,TT,C,nf,Dt,x1,Z1,Z2)
%For numerical solution.

%Z1=2;
%Z2=7e1;
in=1:3:3*C-1;
in2=2:3:3*C-1;
in3=3:3:3*C-1;
%[Prob]= ModelMultiExpo(x,C,TT)
%[Prob]= ModelMultiGammaDeterm(x,C,TT);
[Prob]= ModelMultiGammaDetermDt(x,C,TT,Dt);
Prob(find(~Prob))=realmin;
LogLikelihood=sum(nf.*log(Prob));
%Prior= prod(exp(-(x(in)-1)/Z)/Z)*prod(exp(-(x(in2)-2)/Z)/Z);
%L=-Likelihood*Prior;

LogPrior= sum(log(exp(-(x(in)-x1)/Z1)/Z1))+sum(log(exp(-(x(in2))/Z2)/Z2))-((C-1)*log(1e3));%;
%LogPrior= sum(log(exp(-(x(in)-x1)/Z1)/Z1))+sum(log(exp(-(x(in2))/Z2)/Z2))+sum(log(x(in3)));
%LogPrior=0;
L=-(LogLikelihood+LogPrior);