function [NHuber,PHuber]=Hbr(N,m,s,W,A,theta)
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
            y=Z'*W;
            t=(abs(y)<theta);
            g=y.*t + theta * sign(y).*(1-t);
            dg=sum(t);
            [W,R]=qr(Z*g - W*diag(dg)); %changed here from Z
            W=W*real(inv(W'*W)^(1/2));%Orthonormalization
            crit=abs(sum(W.*Wold)); %chnaged here from W.*
            Niter=Niter+1;
        end
        NHuber=Niter;
        if (NHuber > 999)
              NHuber=-1;
        end
        PHuber=AMARI(W,CC,A);