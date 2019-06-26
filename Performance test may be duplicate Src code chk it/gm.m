function [NGM,PGM]=gm(N,m,s,W,A)
 o=(s'-mean(s',2)*ones(1,N))';%preprocessing part I
epsilon=0.0001;
iter=1000;
Niter=0;
crit=zeros(1,N);
C=cov(o);
CC=C^(-1/2);
Z=CC*o';%Preprocessing part II
W=W*real(inv(W'*W)^(1/2));%FastICA iteration
while(1-min(crit)>epsilon && Niter<iter)
Wold=W;
y=Z'*W;%prewhitened matrix
GM=y.*((1+y.^2).^-2);
dGM=((1+y.^2).^-2);
W=(Z*GM-ones(m,1)*sum(dGM)*W)/N;
W=W*real(inv(W'*W)^(1/2));%Orthonormalization
crit=abs(sum(W.*Wold)); %chnaged here from W.*
Niter=Niter+1;
end
NGM=Niter;
if (NGM > 999)
    NGM=-1;
end
NGM=Niter;
if (NGM > 999)
    NGM=-1;
end

     PGM=AMARI(W,CC,A);
        
