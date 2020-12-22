function [sd,Mean,Skewness, Kurtosis]=extract_his(im)

msize2=size(im);
chunk=zeros(5,5);

sd=0;
Mean=0;
Skewness=0;
Kurtosis=0;
n=1;
for i=1:5:floor(msize2(1)/5)*5
    m=1;
    for j=1:5:floor(msize2(2)/5)*5
    chunk=im(i:i+5-1,j:j+5-1);
    if isempty(find(chunk==0))
        sd(n,m)=std(std(chunk));                        %����Ԫ�صı�׼��
        Mean(n,m)=mean(mean(chunk));                    %����Ԫ�صľ�ֵ  
        Skewness(n,m)=skewness(skewness(chunk));
        Kurtosis(n,m)=kurtosis(kurtosis(chunk));
        m=m+1;
    end
    end
    n=n+1;
end