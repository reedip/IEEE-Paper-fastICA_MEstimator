huber)
y1=wavread('rssd_A.wav');
y2=wavread('rssd_B.wav');
s2=[y1,y2];
[N,m]=size(s2);
%theta=1.345;
theta=1.000;
W=eye(m);
v=s2/chol((s2'*s2)/N);
epsilon=0.0001
iter=1000;
Niter=0;
crit=zeros(1,m);
while(1-min(crit)>epsilon && Niter<iter)
Wold=W;
y=v*W;
t=(abs(y)<theta);
g=y.*t + theta * sign(y).*(1-t);
dg=sum(t);
[W,R]=qr(v'*g - W*diag(dg));
W=W*real(inv(W'*W)^(1/2));
crit=abs(sum(W.*Wold));
Niter=Niter+1;
end



sound(y,16000)
sound(s2,16000)

tanh) N=1
y1=wavread('rssd_A.wav');
y2=wavread('rssd_B.wav');
s2=[y1,y2];
epsilon=0.0001;
iter=1000;
Niter=0;
[N,m]=size(s2);
W=eye(m);
s2=s2/chol((s2'*s2)/N);
C=cov(s2);
CC=C^(-1/2);
Z=CC*s2';
crit=zeros(1,m);
while(1-min(crit)>epsilon && Niter<iter)
Wold=W;
hypTan=tanh(Z'*W);
W=Z*hypTan/N -ones(m,1)*sum(1-hypTan.^2).*W/N;
W=W*real(inv(W'*W)^(1/2));
crit=abs(sum(W.*Wold));
Niter=Niter+1;
end

pow3) N=14
y1=wavread('rssd_A.wav');
y2=wavread('rssd_B.wav');
s2=[y1,y2];
epsilon=0.0001;
iter=1000;
Niter=0;
[N,m]=size(s2);
W=eye(m);
s2=s2/chol((s2'*s2)/N);
C=cov(s2);
CC=C^(-1/2);
Z=CC*s2';
crit=zeros(1,m);
while(1-min(crit)>epsilon && Niter<iter)
Wold=W;
W=(Z*((Z'*W).^3))/N-3*W;
W=W*real(inv(W'*W)^(1/2));
crit=abs(sum(W.*Wold));
Niter=Niter+1;
end

Modified Huber with Preprocessing)N=11
y1=wavread('rssd_A.wav');
y2=wavread('rssd_B.wav');
s2=[y1,y2];
[N,m]=size(s2);
theta=1.345;
%theta=1.000;
W=eye(m);
v=s2/chol((s2'*s2)/N);
epsilon=0.0001
iter=1000;
Niter=0;
crit=zeros(1,m);
C=cov(v);
CC=C^(-1/2);
Z=CC*v';
while(1-min(crit)>epsilon && Niter<iter)
Wold=W;
y=v*W;
t=(abs(y)<theta);
g=y.*t + theta * sign(y).*(1-t);
dg=sum(t);
[W,R]=qr(v'*g - W*diag(dg));
W=W*real(inv(W'*W)^(1/2));
crit=abs(sum(W.*Wold));
Niter=Niter+1;
end

