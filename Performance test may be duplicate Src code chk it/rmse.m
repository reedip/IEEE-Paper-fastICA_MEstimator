function[RMSE]=rmse(estimated,actual)
%The function returns a 1-D root mean square error for the actual and
%estimated signal
%Since ICA's estimated signals can be in a different direction than actual
%signal, therefore RMSE would be for the absolute value of the signals
[N,m]=size(actual);%assuming both estimated and actual signals are of the same size. This would be because we are not worrying about Eigenvalue and eigenvectors
RMSE=zeros(m,1);
actual=actual/norm(actual);
estimated=estimated/norm(estimated);
for j=1:m
    rms=0;
    for i=1:N
        rms=rms+((abs(estimated(i,j)) - abs(actual(i,j))).^2);
    end
    rms=rms/N;
    RMSE(j)=sqrt(rms);
end
%IMP:Normalized the data using the W/norm(W) method.Hope it works
