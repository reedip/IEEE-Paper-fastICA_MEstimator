function [NTanh,PTanh]=htan(N,m,s,W,A)
o=(s'-mean(s',2)*ones(1,N))';
epsilon=0.0001;
iter=1000;
Niter=0;
crit=zeros(1,N);
C=cov(o);
CC=C^(-1/2);
Z=CC*o';
W=W*real(inv(W'*W)^(1/2));%FastICA iteration
while(1-min(crit)>epsilon && Niter<iter)
Wold=W;
hypTan=tanh(Z'*W);
W=Z*hypTan/N -ones(m,1)*sum(1-hypTan.^2).*W/N;
W=W*real(inv(W'*W)^(1/2));
crit=abs(sum(W.*Wold));
Niter=Niter+1;
end
NTanh=Niter;
if (NTanh > 999)
    NTanh=-1;
end
PTanh=AMARI(W,CC,A);