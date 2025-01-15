function alp=alphafind()
global K N  H  P_c P_max zeta w  B  Htill
ep=10^-6;
nvar=2*K;
Aeq=[ones(1,K),zeros(1,K)];  %Aeq*x<=Beq
beq=N;
Ain=[zeros(1,K),ones(1,K)];
bin=P_max;
lb=[ones(K,1);ep*ones(K,1)];
ub=[(N-K)*ones(K,1);P_max*ones(K,1)];
Htill=mean(H')';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%genetic
intvar=[1:K];
%options = gaoptimset('PenaltyFactor',30,'InitialPenalty',100,'StallGenLimit',10,'TolFun',1e-7,...
%             'Generations',200,'SelectionFcn',{@selectiontournament,4},'PopInitRange',[0;P_max],'Display','off',...
%             'PopulationSize',20,'TimeLimit',4);
A=[Ain;Aeq;-Aeq];  % equality constraint can not use in ga, thus model in this way: Ax<b and Ax>b.
b=[bin;beq;-beq];
[xga,etabar,exitflagalpha,outputalpha] = ga(@objalpha,nvar,A,b,[],[],lb,ub,@conalpha,intvar,[]);
NN=xga(1:K)';
P=xga(K+1:2*K)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%fmincon
%options = optimoptions('fmincon','Algorithm','interior-point','GradConstr','on','GradObj','on','Hessian','bfgs');
% X0=[1/K*ones(K,1);P_max/K*ones(K,1)];
% [xf,etabar,exitflagalpha,outputalpha] = fmincon(@objalpha,X0,Ain,bin,Aeq,beq,lb,ub,@conalpha,[]);x=xf';
% NN=xf(1:K);
% NN=round(NN);
% P=xf(K+1:2*K);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%alpha
etabar=-etabar;
B*w.*NN.*log2(1+Htill.*P./NN);
alp=B*w.*NN.*log2(1+Htill.*P./NN)/(etabar*P_c)-zeta*P/P_c;
while abs(sum(alp)-1)>10^-8 || min(alp)<0
    [xga,etabar,exitflagalpha,outputalpha] = ga(@objalpha,nvar,A,b,[],[],lb,ub,@conalpha,intvar,[]);
    NN=xga(1:K)';
    P=xga(K+1:2*K)';
    etabar=-etabar;
    alp=B*w.*NN.*log2(1+Htill.*P./NN)/(etabar*P_c)-zeta*P/P_c;
end

