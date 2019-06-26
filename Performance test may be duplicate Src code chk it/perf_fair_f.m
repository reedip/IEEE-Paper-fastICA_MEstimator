function [N_output,A_output,F_output]=perf_fair_f()
N_output=zeros(100,5);
A_output=zeros(100,5);
F_output=zeros(100,5);
c=0.0;
t=0;
while(c<2.0)
    nsamples=100;
    c=c+0.02;
    t=t+1;
    fprintf('Calculating for c=%d',c); 
while(nsamples<501)
Nf=0;
Af=0;
f2cf=0;
    s1=binary(nsamples);
    s2=root3(nsamples);
    s3=signal(nsamples);
    s=[s1;s2;s3];
    fprintf('Calculating %d samples\n',nsamples);
 for SNRlevel = -20:20
        fprintf('Calculating %d Noise\n',SNRlevel+21);
    for i = 1:100
                fprintf('Calculating %d Iteration\n',i);
               w=eye(3);%10 random weights
        A = 2*rand(3)-1;
        x = (A*s)';
        n = randn(nsamples,3);
        varn = sumsqr(n);
        vary = sumsqr(x);
        n = n/varn*vary/10^(SNRlevel/10);
        x = x + n;
        [N,m]=size(x);
                [N9,A9]=fair(N,m,x,w,A,c); 
        if (N9 ~= -1)
            Nf=Nf+N9;
            Af=Af+A9;
        else
            f2cf=f2cf+1;
        end
   end
 end
      if(f2cf ~= 4100)
        N_output(t,(nsamples/100))=(Nf/(4100 - f2cf));
        A_output(t,(nsamples/100))=(Af/(4100 - f2cf));
        F_output(t,(nsamples/100))=f2cf;
       end 
       csvwrite('N_output.csv',N_output);
       csvwrite('A_output.csv',A_output);
       csvwrite('F_output.csv',F_output);
        nsamples=nsamples+100;
end
end