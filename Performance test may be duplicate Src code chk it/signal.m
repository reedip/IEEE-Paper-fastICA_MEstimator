function [ sqrtN ]= signal(N)
sqrtN=ones(1,N);
for i=1:N
sqrtN(i)=sin(i);
end
%i=1:1000;
%plot(i,sqrtN)
%axis ( [ 0 100 -1 1 ])