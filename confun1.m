function [c,ceq,Gc,Gceq] = confun(x)
%c_i(x)<=0, i=1,...,K,K+1
%ceq_i(x)=0, i=1,...,K
global N R  rho k B nvars H
p(k,:)=zeros(1,N);
nz=logical(rho(k,:));
p(k,nz)=x;
r=find_r(p(k,:));
c=-((B*(sum(rho(k,:).*r(k,:))))-R(k));
ceq=[];

%%%%%%%%%%%%%%%%%%%%%gradian
Gceq=[];
Gc=zeros(nvars,1);
j=1;
for n=1:N
    if rho(k,n)~=0
        Gc(j,1)=-B*rho(k,n)*H(k,n)/((1+H(k,n)*p(k,n))*log(2));
        j=j+1;
    end
end
