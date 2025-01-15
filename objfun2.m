function [f,Gf]=objfun2(x)
 
 
global K N P_c zeta N0 H w k rhu B alpha f2 nvars
%%%%
j=1;
ppp(k,:)=zeros(1,N);
nz=logical(rhu(k,:));
ppp(k,nz)=x;
r=find_r(ppp(k,:));
f=B*(w(k)*sum((rhu(k,:).*r(k,:)))/(alpha(k)*P_c+zeta*sum(ppp(k,:))));
f2=f*(alpha(k)*P_c+zeta*sum(ppp(k,:)));
f=-f;% beacause of the problem is max instead of min
f2=-f2;
%%%%%%%%%%%%%%%%%%%%%%%gradian
Gf=zeros(nvars,1);
j=1;
for n=1:N
    if rhu(k,n)~=0
        Gf(j,1)=(B*w(k)*rhu(k,n)*H(k,n)/((1+H(k,n)*ppp(k,n))*log(2))*...
            (alpha(k)*P_c+zeta*sum(ppp(k,:)))-zeta*B*(w(k)*sum((rhu(k,:).*r(k,:)))))...
            /(alpha(k)*P_c+zeta*sum(ppp(k,:)))^2;
        j=j+1;
    end
end
Gf=-Gf;

