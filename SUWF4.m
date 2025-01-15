function[miu]=SUWF4(k)
global W H K N counter4 counter2 counter3 counterematris counterev R v
miu_int=max(1./H(k,:));
miu(1,k)=miu_int+.7;
%Matris_miu(k,1)=miu(1,k);
miu;
counter4=1;
counter2=1;
counter3=1;
counterematris=0;
counterev=0;

v=[];

