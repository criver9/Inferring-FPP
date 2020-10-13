function [XMIN, ChiMIN, I, xT, ChiT]=FittingHierarMultiGammaModelDeterm(UTimes,c,xprev,xprev2,nf,Dt,TauMin,Z1,Z2)
  N=10;
  xo=10*rand(N,3*c-1);
  in=1:3:3*(c-1)+2;
  in2=2:3:3*(c-1)+2;
  in3=3:3:3*(c-1)+2;
  LB=zeros(1,3*c-1);
  LB(in)=TauMin;
  LB(in2)=0;
  LB(in3)=0;
  
  
  
  
  %xo(:,in)=xo(:,in)*1e3+TauMin;
  %xo(:,in2)=xo(:,in2)*4e3;
  xo(:,in)=xo(:,in)*1+TauMin;
  xo(:,in2)=xo(:,in2)*1;
  
  %LB(in)=1;
  %LB(in2)=2;
  LB
 %options= optimoptions('fmincon','Algorithm','interior-point');
 options = optimset('MaxIter',20000,'MaxFunEvals',20000); 
 for j=1:N
   [xT(j,:),ChiT(j)]=fminsearchbnd(@(x)LogObjectFuncMultiGammaModelDeterm(x,UTimes,c,nf,Dt,TauMin,Z1,Z2),xo(j,:),LB,[],options);
 end
 
 if c==1
  %xinit=[xprev rand+1 rand+2]
  xinit=[xprev rand+TauMin rand]
  xinit2=[xprev rand+TauMin rand]
 else
  %xinit=[xprev 0 1 rand+2]
  xinit=[xprev 0 rand+TauMin rand];
  xinit2=[xprev 1e-4 rand+TauMin rand];
  
  
 end 
  r=rand(5,3*c-1)*1e-4
  [xT(N+1,:) ChiT(N+1)]=fminsearchbnd(@(x)LogObjectFuncMultiGammaModelDeterm(x,UTimes,c,nf,Dt,TauMin,Z1,Z2),xinit,LB,[],options);
  [xT(N+2,:) ChiT(N+2)]=fminsearchbnd(@(x)LogObjectFuncMultiGammaModelDeterm(x,UTimes,c,nf,Dt,TauMin,Z1,Z2),xinit2,LB,[],options);
  [xT(N+3,:) ChiT(N+3)]=fminsearchbnd(@(x)LogObjectFuncMultiGammaModelDeterm(x,UTimes,c,nf,Dt,TauMin,Z1,Z2),xprev2{1,c},LB,[],options);
  [xT(N+4,:) ChiT(N+4)]=fminsearchbnd(@(x)LogObjectFuncMultiGammaModelDeterm(x,UTimes,c,nf,Dt,TauMin,Z1,Z2),xprev2{1,c}+r(1,:),LB,[],options);
  
  ChiT
  [ChiMin I]=min(ChiT)
  XMin=xT(I,:);
  
  [XMIN, ChiMIN]=fminsearchbnd(@(x)LogObjectFuncMultiGammaModelDeterm(x,UTimes,c,nf,Dt,TauMin,Z1,Z2),XMin,LB,[],options);
  
  