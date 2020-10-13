% Estimating Hessians.

s=1;
DiagMatrix=cell(1,5);
for m=2:18
k=(m-2)*2+1;
Xminm=XminCanonDim{m-1,1};
%Xminm=XminCanonN;
[hess,err]=hessian(@(x)MinusLogLikelihood(x,Spiked,TFuture,Sigma,m,uPrior),Xminm)
[Vm,Dm]=eig(hess);
DiagMatrix{1,s}=Dm;
Det(m-1)=prod(diag(Dm));
Evidence(s)=-ChiMinCanon(m-1,1)*(2*pi)^(k/2)*sqrt((1/prod(diag(Dm))));
clear  Xminm Vm Dm
s=s+1;
end
%FIRST=log(-ChiMinCanon(:,2));
%K=1:2:25;
%SECOND=K/2*log(2*pi);
%Det(7)=[];
%THIRD=1/2*log(real(Det));





