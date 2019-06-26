function [NCau,PCau]=cau(N,m,s,W,A,ck)
    o=(s'-mean(s',2)*ones(1,N))';%preprocessing part I
    epsilon=0.0001;
    iter=1000;
    Niter=0;
    crit=zeros(1,N);
    C=cov(o);
    CC=C^(-1/2);
    Z=CC*o';%Preprocessing part II
    %ck=0.225+(2.05-0.225)*rand;%Randomly Generated Cauchy COnstant
    %ck=0.0594+(0.1706-0.0594)*rand;
    W=W*real(inv(W'*W)^(1/2));%FastICA iteration
    while(1-min(crit)>epsilon && Niter<iter)
        Wold=W;
        y=(Z'*W)/ck;%prewhitened matrix
        dCauchy=((1+y.^2).^-2);
        Cauchy=y.*dCauchy;
        W=(Z*Cauchy-ones(m,1)*sum(dCauchy)*W)/N;
        W=W*real(inv(W'*W)^(1/2));%Orthonormalization
        crit=abs(sum(W.*Wold)); %chnaged here from W.*
        Niter=Niter+1;
    end
     NCau=Niter;
    if (NCau > 999)
         NCau=-1;
    end
     PCau=AMARI(W,CC,A);
        