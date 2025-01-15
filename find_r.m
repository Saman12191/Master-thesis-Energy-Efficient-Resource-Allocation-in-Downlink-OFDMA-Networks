function r=find_r(p)
global W N0 H N K k
r(k,:)=zeros(1,N);

    for n=1:N
        r(k,n)=log2(1+(p(n)*(H(k,n))));
    end

        