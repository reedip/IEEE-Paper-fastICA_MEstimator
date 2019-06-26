function [ sqrtN ]= root3(N)
sqrtN=ones(1,N).*sqrt(3);
a=-1;
for i=2:N
sqrtN(i)=sqrtN(i-1).*a;
end
%i=1:1000;
%plot(i,sqrtN)
%axis ( [ 0 100 -sqrt(3) sqrt(3) ])