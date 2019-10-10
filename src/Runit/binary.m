function [ sqrtN ]= binary(N)
sqrtN=ones(1,N);
a=-1;
for i=2:N
sqrtN(i)=sqrtN(i-1).*a;
end
%i=1:1000;
%plot(i,sqrtN)
%axis ( [ 0 100 -1 1 ])