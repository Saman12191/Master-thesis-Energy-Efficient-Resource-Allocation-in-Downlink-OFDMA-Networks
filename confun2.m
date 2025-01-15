function [c,ceq,Gc,Gceq] = confun2(x)
%c_i(x)<=0, i=1,...,K,K+1
%ceq_i(x)=0, i=1,...,K
global N R  rhu k B nvars H
ppp(k,:)=zeros(1,N);
nz=logical(rhu(k,:));
ppp(k,nz)=x;
r=find_r(ppp(k,:));
c=-((B*(sum(rhu(k,:).*r(k,:))))-R(k));
ceq=[];

%%%%%%%%%%%%%%%%%%%%%gradian
Gceq=[];
Gc=zeros(nvars,1);
j=1;
for n=1:N
    if rhu(k,n)~=0
        Gc(j,1)=-B*rhu(k,n)*H(k,n)/((1+H(k,n)*ppp(k,n))*log(2));
        j=j+1;
    end
end
