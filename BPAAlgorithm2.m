function[eta_EE,P_1,iter]=BPAAlgorithm2(k)
global ro B alpha R w  zeta K N P_max P_c e H
S=ro;

Matris_miuu=[];
Matris_Ptot=[];
Matris_v=[];
count4=1;
alpha;
iter=0;
for count4=1
    [miu]=SUWF4(k);
    [miu,phat]=SUWF3(miu,k,S);
    cap=0;
    P_1=phat;
    P1=sum((P_1));
    cap=w(k,1)*R(k,1);
    d=w(k,1)*B*log2(e)/(miu(1,k));
    eta_EE_1=cap/(zeta*P1+alpha(k,1)*P_c);
    eta_SE_1=cap;
    d_1=d;
    if eta_EE_1>=d_1/zeta
        P=P_1;
        eta_EE=eta_EE_1;
        eta_SE=eta_SE_1;
        break
        disp('return inja bood')
    else
        P_2=P_1;
        P2=P1;
        kk=1.5;
        P1=min(kk*P1,P_max); %jame phat ha
        P_inp=P1;
        [P_1,miuu]=MLWF(P1,phat,H,S,k);
        P=sum(P_1);
        P_out=P;%%%jame P*ha
        [d_1,eta_EE_1,eta_SE_1]=compare(S,P_1,miuu,k,P);
    end
    %aa=1;
    t=1;
    while   eta_EE_1<d_1/zeta & P<P_max
        t=t+1;
        P_2=P_1;
        P2=P1;
        P1=min(kk*P1,P_max);
        P_inp=P1;
        [P_1,miuu]=MLWF(P1,phat,H,S,k);
        P=sum(P_1);
        P_out=P;%%%jame P*ha
        P_diff=P_out-P_inp;
        P=sum(P_1);
        [d_1,eta_EE_1,eta_SE_1]=compare(S,P_1,miuu,k,P);
    end
    
    if eta_EE_1<=d_1/zeta
        P=P1;
        eta_EE=eta_EE_1;
        eta_SE=eta_SE_1;
        return
    end
    tt=t;
    while abs(eta_EE_1-d_1/zeta)>eps
        if tt>10000
            break
        end
        val=abs(eta_EE_1-d_1/zeta);
        P2=sum(P_2);
        P=(P1+P2)/2;
        P1=P;
        P_inp=P1;
        [P_1,miuu]=MLWF(P1,phat,H,S,k);
        P=sum(P_1);
        P_out=P;%%%jame P*ha
        P_diff=P_out-P_inp;
        [d_1,eta_EE_1,eta_SE_1]=compare(S,P_1,miuu,k,P);
        if eta_EE_1<d_1/zeta
            P_2=P_1; P2=P;
        else
            P_1=P_1; P1=P;
        end
    tt=tt+1;
    end
    iter=tt;
    eta_EE=eta_EE_1;
    eta_SE=eta_SE_1;
    P_1=P_1;
end

