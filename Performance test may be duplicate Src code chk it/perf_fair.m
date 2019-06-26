%M-Estimator-Evaluation-Function%
function [N_output,A_output,F_output]=perf_fair()
nsamples=500;
N_output=zeros(10,9);
A_output=zeros(10,9);
F_output=zeros(10,9);
while(nsamples<5001)
    empty=0;
Np=0;
Ap=0;
f2cp=0;
Nc1=0;
Ac1=0;
f2cc1=0;
Nw1=0;
Aw1=0;
f2cw1=0;
Ng=0;
Ag=0;
f2cg=0;
Nt=0;
At=0;
f2ct=0;
Nh=0;
Ah=0;
f2ch=0;
Nc2=0;
Ac2=0;
f2cc2=0;
Nw2=0;
Aw2=0;
f2cw2=0;
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
        [N1,A1]=Hbr(N,m,x,w,A);
        if (N1 ~= -1)
                Nh=Nh+N1;
                Ah=Ah+A1;
        else
			f2ch=f2ch+1;
        end
        [N2,A2]=cau(N,m,x,w,A);
        if (N2 ~= -1)
			Nc1=Nc1+N2;
			Ac1=Ac1+A2;
		else
			f2cc1=f2cc1+1;
        end
        [N3,A3]=wel(N,m,x,w,A);
		if (N3 ~= -1)
			Nw1=Nw1+N3;
			Aw1=Aw1+A3;
		else
			f2cw1=f2cw1+1;
		end
       	[N4,A4]=gm(N,m,x,w,A);
		if (N4 ~= -1)
			Ng=Ng+N4;
			Ag=Ag+A4;
		else
			f2cg=f2cg+1;
		end
        [N5,A5]=htan(N,m,x,w,A);
		if (N5 ~= -1)
			Nt=Nt+N5;
			At=At+A5;
		else
			f2ct=f2ct+1;
		end
        [N6,A6]=pow3(N,m,x,w,A);
		if (N6 ~= -1)
			Np=Np+N6;
			Ap=Ap+A6;
		else
			f2cp=f2cp+1;
		end
       		[N7,A7]=cau1(N,m,x,w,A);
		if (N7 ~= -1)
			Nc2=Nc2+N7;
			Ac2=Ac2+A7;
		else
			f2cc2=f2cc2+1;
		end
        [N8,A8]=wel1(N,m,x,w,A);
		if (N8 ~= -1)
			Nw2=Nw2+N8;
			Aw2=Aw2+A8;
		else
			f2cw2=f2cw2+1;
		end
        		[N9,A9]=fair(N,m,x,w,A);  
		if (N9 ~= -1)
        	Nf=Nf+N9;
			Af=Af+A9;
            empty=empty+1;
		else
			f2cf=f2cf+1;
        end
   end
 end
  empty
       if (f2ch ~= 4100)
        N_output((nsamples/500),1)=(Nh/(4100- f2ch));
        A_output((nsamples/500),1)=(Ah/(4100- f2ch));
        F_output((nsamples/500),1)=f2ch;
       end 
      if(f2cc1 ~= 4100)
        N_output((nsamples/500),2)=(Nc1/(4100 - f2cc1));
        A_output((nsamples/500),2)=(Ac1/(4100 - f2cc1));
        F_output((nsamples/500),2)=f2cc1;
       end 
      if(f2cw1 ~= 4100)
        N_output((nsamples/500),3)=(Nw1/(4100 - f2cw1));
        A_output((nsamples/500),3)=(Aw1/(4100 - f2cw1));
        F_output((nsamples/500),3)=f2cw1;
       end 
      if(f2cg ~= 4100)
        N_output((nsamples/500),4)=(Ng/(4100 - f2cg));
        A_output((nsamples/500),4)=(Ag/(4100 - f2cg));
        F_output((nsamples/500),4)=f2cg;
       end 
      if(f2ct ~= 4100)
        N_output((nsamples/500),5)=(Nt/(4100 - f2ct));
        A_output((nsamples/500),5)=(At/(4100 - f2ct));
        F_output((nsamples/500),5)=f2ct;
       end 
      if(f2cp ~= 4100)
        N_output((nsamples/500),6)=(Np/(4100 - f2cp));
        A_output((nsamples/500),6)=(Ap/(4100 - f2cp));
        F_output((nsamples/500),6)=f2cp;
       end 
      if(f2cc2 ~= 4100)
        N_output((nsamples/500),7)=(Nc2/(4100 - f2cc2));
        A_output((nsamples/500),7)=(Ac2/(4100 - f2cc2));
        F_output((nsamples/500),7)=f2cc2;
       end
      if(f2cw2 ~= 4100)
        N_output((nsamples/500),8)=(Nw2/(4100 - f2cw2));
        A_output((nsamples/500),8)=(Aw2/(4100 - f2cw2));
        F_output((nsamples/500),8)=f2cw2;
       end 
      if(f2cf ~= 4100)
        N_output((nsamples/500),9)=(Nf/(4100 - f2cf));
        A_output((nsamples/500),9)=(Af/(4100 - f2cf));
        F_output((nsamples/500),9)=f2cf;
       end  
       csvwrite('N_output.csv',N_output);
       csvwrite('A_output.csv',A_output);
       csvwrite('F_output.csv',F_output);
        nsamples=nsamples+500;
end
