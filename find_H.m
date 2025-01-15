function h=find_H()
global N0 K N W B
nTap = 6;
h=zeros(K,N);
l=0;

for k=1:K
    ht =1/sqrt(nTap)*(randn(nTap,1) + j*randn(nTap,1));
    h(k,:)=(abs(fft(ht,N)).').^2/N0/B;
end
