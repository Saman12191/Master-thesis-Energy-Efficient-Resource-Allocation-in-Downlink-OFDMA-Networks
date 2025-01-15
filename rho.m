function xx=testrofix()
clear all
clc
global K N N0 H W P_c P_max zeta w R S
K=4;
N=2;
N0=1e-6;
zeta=.38;
P_c=10;
P_max=20;
R=100000*ones(K,1);
W=1.8*10^6;
w=ones(K,1);
H=find_H();
nvars=K*N;%Number of variables
lb=zeros(K*N,1);%lower bound for x
ub=inf*ones(K*N,1);%upper bound for x
%%% linear inequality constraint: A.x<=b
A=ones(1,K*N);
b=P_max;
e=3.1416;
B=1.08*1e6;
alpha=.25*ones(1,4);
eta_EE=zeros(1,K);
%j=1;
%for P_max=2:1:8
for j=0:1:5
    j
    roTotal=ones(1,72);
    ro=zeros(K,N);
    g=h+20*j;
    g_mean(1,j+1)=floor(mean2(g));
    for k=1:K
        [val,col]=min(roTotal.*g(k,:));
        nhat=col;
        roTotal(1,nhat)=1000;
        ro(k,nhat)=1;
    end
    for kk=1:72
        if roTotal(1,kk)==1000
            roTotal(1,kk)=0;
        end
    end
    S=ro;
    xx=ga(@objfunro,nvars,A,b,[],[],lb,ub,@confunro)
    cap=zeros(1,K);
    for k=1:K
        P=0;
        for n=1:N
            if ro(k,n)>0
                cap(1,k)=cap(1,k)+log2(1+xx(1,n)*g(k,n));  %% khoruji algorithm zhenetic
                P=P+xx(1,n);
            end
        end
        eta_EE(1,k)=W*w(1,k)*cap(1,k)/(zeta*P+alpha(1,k)*P_c);
        eta_SE(1,k)=W*w(1,k)*cap(1,k);
        
    end
    l=find(roTotal==0);
    m=length(l);
    disp('19')
    while m<72
        m;
        q=min(eta_EE);
        [val,col]=min(eta_EE);
        khat=col;
        [val,col]=max(roTotal.*g(khat,:));
        ntild=col;
        roTotal(1,ntild)=0;
        ro(khat,ntild)=1;
        S=ro;
        xx=feval(@pso2,S,g,P_max);
        cap=zeros(1,K);
        for k=1:K
            P=0;
            for n=1:N
                if ro(k,n)>0
                    cap(1,k)=cap(1,k)+log2(1+xx(1,n)*g(k,n));
                    P=P+xx(1,n);
                end
            end
            
            eta_EE(1,k)=W*w(1,k)*cap(1,k)/(zeta*P+alpha(1,k)*P_c);
            eta_SE(1,k)=W*w(1,k)*cap(1,k);
            
        end
        l=find(roTotal==0);
        m=length(l);
    end
    eta_EEE(1,j+1)=mean(eta_EE);
    eta_SEE(1,j+1)=mean(eta_SE);
    %Pow(1,j+1)=P_max;
    %j=j+1;
end
%plot(Pow,eta_EEE,'-*','linewidth',2,'markersize',8)
% eta_EEE(1,j+1)=mean(eta_EE);
plot(floor(10*log10(g_mean)),eta_EEE,'-*','linewidth',2,'markersize',8)
figure
plot(floor(10*log10(g_mean)),eta_SEE,'-*','linewidth',2,'markersize',8)
profile report
end
%%%%%
function [c,ceq] = confunro(x)
%c_i(x)<=0, i=1,...,K,K+1
%ceq_i(x)=0, i=1,...,K
global K N R P_max S
rho=S;
p=reshape(x,K,N);
r=find_r(p);
for i=1:K 
    c(i)=-(sum(rho(i,:).*r(i,:))-R(i));
end
% other constraints are linear and we wright these as forme A.x<=b and
% Aeq.x=beq in test file
%    c(K+1)=P_max-sum(sum(p));
%for i=1:K
%    ceq(i)=sum(rho(:,i))-1;
%end
ceq=[];
end
function f=objfunro(x)
%x=[rho_(1,1),...,rho_(1,N),rho_(2,1),...,rho_(2,N),...,rho_(K,1),...,rho_(K,N),
%p_(1,1),...,p_(1,N),p_(2,1),...,p_(2,N),...,p_(K,1),...,p_(K,N)]. Note: x
%is a vertical vector with 2*K*N element. 
global K N P_c zeta N0 H w S
%%%%
rho=S;
p=reshape(x,K,N);
%%%%
r=find_r(p);
f=sum(sum(diag(w)*(rho.*r)))/(P_c+zeta*sum(sum(p)));
f=-f;% beacause of the problem is max instead of min
end
