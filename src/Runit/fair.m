function [Nfair,Pfair]=fair(N,m,s,W,A,c)
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
            dg=(1+(abs(y)/c));
            g=y./(dg);
            %W=(Z*g'-ones(m,1)*sum(dg)*W)/N;
            [W,R]=qr(Z*g - ones(m,1)*sum(dg)*W);
            W=W*real(inv(W'*W)^(1/2));%Orthonormalization
            crit=abs(sum(W.*Wold)); %chnaged here from W.*
            Niter=Niter+1;
            %fprintf('%d:',Niter);
        end
        Nfair=Niter;
        if (Nfair > 999)
                Nfair=-1;
        end
        Pfair=AMARI(W,CC,A);
        
