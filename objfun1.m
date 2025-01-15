function [f,Gf]=objfun(x)
 
 
global K N P_c zeta N0 H w k rho B alpha f2 nvars
%%%%
j=1;
p(k,:)=zeros(1,N);
nz=logical(rho(k,:));
p(k,nz)=x;
r=find_r(p(k,:));
f=B*(w(k)*sum((rho(k,:).*r(k,:)))/(alpha(k)*P_c+zeta*sum(p(k,:))));
f2=f*(alpha(k)*P_c+zeta*sum(p(k,:)));
f=-f;% beacause of the problem is max instead of min
f2=-f2;
%%%%%%%%%%%%%%%%%%%%%%%gradian
Gf=zeros(nvars,1);
j=1;
for n=1:N
    if rho(k,n)~=0
        Gf(j,1)=(B*w(k)*rho(k,n)*H(k,n)/((1+H(k,n)*p(k,n))*log(2))*...
            (alpha(k)*P_c+zeta*sum(p(k,:)))-zeta*B*(w(k)*sum((rho(k,:).*r(k,:)))))...
            /(alpha(k)*P_c+zeta*sum(p(k,:)))^2;
        j=j+1;
    end
end
Gf=-Gf;

