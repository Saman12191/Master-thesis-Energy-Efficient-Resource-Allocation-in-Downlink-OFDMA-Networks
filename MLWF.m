function [P_1,miuu]=MLWF(P1,phat,g,S,k)
global N K w 
miuu=0;
sum1=0; 
sum2=0;
sum3=0;
b=zeros(1,N);
c=zeros(1,N);
pstar=zeros(1,N);
    for n=1:N
        if S(k,n)>0
         sum3=sum3+w(k,1);
         sum1=sum1+phat(1,n);
         sum2=sum2+1/g(k,n)+phat(1,n);
        end
    end
miuu=(P1-sum1+sum2)/sum3;
    for n=1:N
        if S(k,n)>0
        c(1,n)=(miuu-(1/g(k,n)+phat(1,n)));
        b(1,n)=max(c(1,n),0);
        pstar(1,n)=phat(1,n)+b(1,n);
        end
    end
P_1=pstar;


