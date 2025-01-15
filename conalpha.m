function [cc,cceq]=conalpha(x)
global K  B  R Htill
NN=x(1:K);
P=x(K+1:2*K);
cc=zeros(1,K);
for k=1:K
    cc(k)=R(k)-NN(k)*B*log2(1+Htill(k)*P(k)/NN(k));
end
cceq=[];
end