function [huber,cau_1,cau2,tanh,wel_1,wel2,pow,fair1,gm1]=perf_fair()
N1=zeros(100,9);
N2=zeros(100,9);
N3=zeros(100,9);
N4=zeros(100,9);
N5=zeros(100,9);
N6=zeros(100,9);
N7=zeros(100,9);
N8=zeros(100,9);
N9=zeros(100,9);
A1=zeros(100,1);
A2=zeros(100,1);
A3=zeros(100,1);
A4=zeros(100,1);
A5=zeros(100,1);
A6=zeros(100,1);
A7=zeros(100,1);
A8=zeros(100,1);
A9=zeros(100,1);
a_mean1=zeros(61,9);
a_mean2=zeros(61,9);
a_mean3=zeros(61,9);
a_mean4=zeros(61,9);
a_mean5=zeros(61,9);
huber=zeros(61,5);
cau_1=zeros(61,5);
cau2=zeros(61,5);
tanh=zeros(61,5);
wel_1=zeros(61,5);
wel2=zeros(61,5);
pow=zeros(61,5);
fair1=zeros(61,5);
gm1=zeros(61,5);
nsamples=5000;
while(nsamples<5001)
    s1=binary(nsamples);
    s2=root3(nsamples);
    s3=signal(nsamples);
    s=[s1;s2;s3];
    fprintf('Calculating %d samples\n',nsamples);
 for SNRlevel = -30:30;
        fprintf('Calculating %d Noise\n',SNRlevel+31);
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
        [N1(i,1),A1(i,1)]=Hbr(N,m,x,w,A);
        [N2(i,2),A2(i,1)]=cau(N,m,x,w,A);
        [N3(i,3),A3(i,1)]=wel(N,m,x,w,A);
        [N4(i,4),A4(i,1)]=gm(N,m,x,w,A);
        [N5(i,5),A5(i,1)]=htan(N,m,x,w,A);
        [N6(i,6),A6(i,1)]=pow3(N,m,x,w,A);
        [N7(i,7),A7(i,1)]=cau1(N,m,x,w,A);
        [N8(i,8),A8(i,1)]=wel1(N,m,x,w,A);
        [N9(i,9),A9(i,1)]=fair(N,m,x,w,A);  
   end
        if(nsamples<1001)
        a_mean1(SNRlevel+31,1) = mean(A1);
        a_mean1(SNRlevel+31,2) = mean(A2);
        a_mean1(SNRlevel+31,3) = mean(A3);
        a_mean1(SNRlevel+31,4) = mean(A4);
        a_mean1(SNRlevel+31,5) = mean(A5);
        a_mean1(SNRlevel+31,6) = mean(A6);
        a_mean1(SNRlevel+31,7) = mean(A7);
        a_mean1(SNRlevel+31,8) = mean(A8);
        a_mean1(SNRlevel+31,9) = mean(A9);
        end
        if((nsamples>1001)&&(nsamples<2001))
        a_mean2(SNRlevel+31,1) = mean(A1);
        a_mean2(SNRlevel+31,2) = mean(A2);
        a_mean2(SNRlevel+31,3) = mean(A3);
        a_mean2(SNRlevel+31,4) = mean(A4);
        a_mean2(SNRlevel+31,5) = mean(A5);
        a_mean2(SNRlevel+31,6) = mean(A6);
        a_mean2(SNRlevel+31,7) = mean(A7);
        a_mean2(SNRlevel+31,8) = mean(A8);
        a_mean2(SNRlevel+31,9) = mean(A9);
        end
        if((nsamples>2001)&&(nsamples<3001))
        a_mean3(SNRlevel+31,1) = mean(A1);
        a_mean3(SNRlevel+31,2) = mean(A2);
        a_mean3(SNRlevel+31,3) = mean(A3);
        a_mean3(SNRlevel+31,4) = mean(A4);
        a_mean3(SNRlevel+31,5) = mean(A5);
        a_mean3(SNRlevel+31,6) = mean(A6);
        a_mean3(SNRlevel+31,7) = mean(A7);
        a_mean3(SNRlevel+31,8) = mean(A8);
        a_mean3(SNRlevel+31,9) = mean(A9);
        end
        if((nsamples>3001)&&(nsamples<4001))
        a_mean4(SNRlevel+31,1) = mean(A1);
        a_mean4(SNRlevel+31,2) = mean(A2);
        a_mean4(SNRlevel+31,3) = mean(A3);
        a_mean4(SNRlevel+31,4) = mean(A4);
        a_mean4(SNRlevel+31,5) = mean(A5);
        a_mean4(SNRlevel+31,6) = mean(A6);
        a_mean4(SNRlevel+31,7) = mean(A7);
        a_mean4(SNRlevel+31,8) = mean(A8);
        a_mean4(SNRlevel+31,9) = mean(A9);
        end
        if((nsamples>4001)&&(nsamples<5001))
        a_mean5(SNRlevel+31,1) = mean(A1);
        a_mean5(SNRlevel+31,2) = mean(A2);
        a_mean5(SNRlevel+31,3) = mean(A3);
        a_mean5(SNRlevel+31,4) = mean(A4);
        a_mean5(SNRlevel+31,5) = mean(A5);
        a_mean5(SNRlevel+31,6) = mean(A6);
        a_mean5(SNRlevel+31,7) = mean(A7);
        a_mean5(SNRlevel+31,8) = mean(A8);
        a_mean5(SNRlevel+31,9) = mean(A9);
        end
        %plot(10*log10(a_mean));
        huber=[a_mean1(:,1),a_mean2(:,1),a_mean3(:,1),a_mean4(:,1),a_mean5(:,1)];
        cau_1=[a_mean1(:,2),a_mean2(:,2),a_mean3(:,2),a_mean4(:,2),a_mean5(:,2)];
        cau2=[a_mean1(:,3),a_mean2(:,3),a_mean3(:,3),a_mean4(:,3),a_mean5(:,3)];
        tanh=[a_mean1(:,4),a_mean2(:,4),a_mean3(:,4),a_mean4(:,4),a_mean5(:,4)];
        wel_1=[a_mean1(:,5),a_mean2(:,5),a_mean3(:,5),a_mean4(:,5),a_mean5(:,5)];
        wel2=[a_mean1(:,6),a_mean2(:,6),a_mean3(:,6),a_mean4(:,6),a_mean5(:,6)];
        pow=[a_mean1(:,7),a_mean2(:,7),a_mean3(:,7),a_mean4(:,7),a_mean5(:,7)];
        fair1=[a_mean1(:,8),a_mean2(:,8),a_mean3(:,8),a_mean4(:,8),a_mean5(:,8)];
        gm1=[a_mean1(:,9),a_mean2(:,9),a_mean3(:,9),a_mean4(:,9),a_mean5(:,9)];
       csvwrite('huber.csv',huber);
       csvwrite('cau_1.csv',cau_1);
       csvwrite('cau2.csv',cau2);
       csvwrite('wel_1.csv',wel_1);
       csvwrite('wel2.csv',wel2);
       csvwrite('pow3.csv',pow);
       csvwrite('fair1.csv',fair1);
       csvwrite('gm.csv',gm1);
       csvwrite('tanh.csv',tanh);
  end
        nsamples=nsamples+1000;
      
end
 % csvwrite('A1.csv',A1);
 % csvwrite('A2.csv',A2);
 % csvwrite('A3.csv',A3);
 % csvwrite('A4.csv',A4);
 % csvwrite('A5.csv',A5);
 % csvwrite('A6.csv',A6);
 % csvwrite('A7.csv',A7);
 % csvwrite('A8.csv',A8);
 % csvwrite('A9.csv',A9);
  

