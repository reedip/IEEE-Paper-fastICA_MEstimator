function []=perf_fair_f()
HN_output=zeros(100,5);
HA_output=zeros(100,5);
HF_output=zeros(100,5);
CN_output=zeros(100,5);
CA_output=zeros(100,5);
CF_output=zeros(100,5);
WN_output=zeros(100,5);
WA_output=zeros(100,5);
WF_output=zeros(100,5);
c=0.0;
t=0;
while(c<2.0)
    nsamples=100;
    c=c+0.02;
    t=t+1;
    fprintf('Calculating for c=%d',c); 
while(nsamples<501)
Nc=0;
Ac=0;
f2cc=0;
Nh=0;
Ah=0;
f2ch=0;
Nw=0;
Aw=0;
f2cw=0;
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
%         [N1,A1]=Hbr(N,m,x,w,A,c); 
%         if (N1 ~= -1)
%             Nh=Nh+N1;
%             Ah=Ah+A1;
%         else
%             f2ch=f2ch+1;
%         end
         [N3,A3]=fair(N,m,x,w,A,c); 
         if (N3 ~= -1)
             Nw=Nw+N3;
             Aw=Aw+A3;
         else
             f2cw=f2cw+1;
         end
    %	c1=c*3/2;
 %       [N2,A2]=cau(N,m,x,w,A,c1); 
 %       if (N2 ~= -1)
 %           Nc=Nc+N2;
 %           Ac=Ac+A2;
 %       else
 %           f2cc=f2cc+1;
%       end

   end
 end
  %    if(f2ch ~= 4100)
  %      HN_output(t,(nsamples/100))=(Nh/(4100 - f2ch));
  %      HA_output(t,(nsamples/100))=(Ah/(4100 - f2ch));
  %      HF_output(t,(nsamples/100))=f2ch;
  %     end 
  %     if(f2cc ~= 4100)
  %       CN_output(t,(nsamples/100))=(Nc/(4100 - f2cc));
  %       CA_output(t,(nsamples/100))=(Ac/(4100 - f2cc));
  %       CF_output(t,(nsamples/100))=f2cc;
  %     end 
        if(f2cw ~= 4100)
         WN_output(t,(nsamples/100))=(Nw/(4100 - f2cw));
         WA_output(t,(nsamples/100))=(Aw/(4100 - f2cw));
         WF_output(t,(nsamples/100))=f2cw;
        end 
   %    csvwrite('HN_output.csv',HN_output);
   %    csvwrite('HA_output.csv',HA_output);
   %    csvwrite('HF_output.csv',HF_output);
   %    csvwrite('CN_output.csv',CN_output);
   %     csvwrite('CA_output.csv',CA_output);
   %     csvwrite('CF_output.csv',CF_output);
        csvwrite('FN_output.csv',WN_output);
        csvwrite('FA_output.csv',WA_output);
        csvwrite('FF_output.csv',WF_output);
        nsamples=nsamples+100;
end
end