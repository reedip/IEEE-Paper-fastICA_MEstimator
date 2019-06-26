function[NWel,PWel]=wel(N,m,s,W,A)
o=(s'-mean(s',2)*ones(1,N))';%preprocessing part I
epsilon=0.0001;
iter=1000;
Niter=0;
crit=zeros(1,N);
C=cov(o);
CC=C^(-1/2);
Z=CC*o';%Preprocessing part II
%c=2.9846;%Cauchy COnstant
c=1.154+(1.9-1.154)*rand;
%find out welsch constant on the internet. it is required as it may be used
%in several places and the absoulte value is required for proper
%convergence of the algorithm. We can also do an independent study of how
%the Welsch's constant drives FastICA
W=W*real(inv(W'*W)^(1/2));%FastICA iteration
while(1-min(crit)>epsilon && Niter<iter)
Wold=W;
y=Z'*W;%prewhitened matrix
dWelsch=(exp(-(((y)/c).^2)));
Welsch=y.*dWelsch;
W=(Z*Welsch-ones(m,1)*sum(dWelsch)*W)/N;
W=W*real(inv(W'*W)^(1/2));%Orthonormalization
crit=abs(sum(W.*Wold)); %chnaged here from W.*
Niter=Niter+1;
end
NWel=Niter;
if (NWel > 999)
    NWel=-1;
end
     PWel=AMARI(W,CC,A);


