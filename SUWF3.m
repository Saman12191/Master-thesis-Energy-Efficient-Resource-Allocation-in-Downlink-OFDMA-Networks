function[miu,phat]=SUWF3(miu,k,S)
global W H K N counter4 counter2 counter3 counterematris counterev R v
format long
phatt=zeros(1,N);
for n=1:N
    if S(k,n)>0
        matris_max(1,n)=max((miu(1,k)-1/H(k,n)),0);
        phatt(1,n)=max((miu(1,k)-1/H(k,n)),0);
    end
end
phatt;
i=1;
for n=1:N
    if S(k,n)>0
        if phatt(1,n)>0
            v(1,i)=H(k,n);
            counterev=counterev+1;
            Matris_v(counterev)=v(1,i);
            i=i+1;
        end
    end
end
%[f,c,r]=find(v);
%l=length(r);
l=length(v);
t(1,k)=i;
x=(prod(v));
x;
dd=1/l;
counter3=counter3+1;
miu(1,k)=((2^((R(k,1)*N/W)))/x)^dd;
Matris_miu(k,counter3)=miu(1,k);
phat=phatt;
P_tot=sum(phatt);
counterematris=counterematris+1;
Matris_Ptot(k,counterematris)=P_tot;
%Matris_phat(counter2,:)=phat;

counter2=counter2+1;
% if abs(Matris_Ptot(k,counterematris)-.2*P_max)>.01
%     return
% end
while counterematris~=3
    [miu,phat]=SUWF3(miu,k,S);
end

