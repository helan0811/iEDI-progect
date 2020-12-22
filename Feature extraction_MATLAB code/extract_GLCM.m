function [contrast,correlation,energy,homogeneity,entropy]=extract_GLCM(im)

dis=1;
im=double(im);
numLevel=5;

msize2=size(im);
chunk=zeros(5,5);

Contrast=zeros(1,1,4);
Correlation=zeros(1,1,4);
Energy=zeros(1,1,4);
Homogeneity=zeros(1,1,4);
Entropy=zeros(1,1,4);
n=1;
for i=1:5:floor(msize2(1)/5)*5
    m=1;
    for j=1:5:floor(msize2(2)/5)*5
    chunk=im(i:i+5-1,j:j+5-1);
    if isempty(find(chunk==0))
        GLCM2 = graycomatrix(chunk,'GrayLimits',[min(chunk(:)), max(chunk(:))],'NumLevels',numLevel,'Offset',[0 dis; -dis dis; -dis 0; -dis -dis]);
        stats= graycoprops(GLCM2,'all'); 
        for step=1:4
            Contrast(n,m,step)=stats.Contrast(step);
            Correlation(n,m,step)=stats.Correlation(step);
            Energy(n,m,step)=stats.Energy(step);
            Homogeneity(n,m,step)=stats.Homogeneity(step);
            tGLCM=normalizeGLCM(GLCM2(:,:,step));
            stats.Entropy(step)=calculateEntropy(tGLCM);
            Entropy(n,m,step)=stats.Entropy(step);
        end
        
        m=m+1;
    end
    end
    n=n+1;
end


[n1,m1,k1]=size(Contrast);
contrast=zeros(n1,m1);
correlation=zeros(n1,m1);
energy=zeros(n1,m1);
homogeneity=zeros(n1,m1);
entropy=zeros(n1,m1);
contrast(:,:)=(Contrast(:,:,1)+Contrast(:,:,2)+Contrast(:,:,3)+Contrast(:,:,4))/4;
correlation(:,:)=(Correlation(:,:,1)+Correlation(:,:,2)+Correlation(:,:,3)+Correlation(:,:,4))/4;
energy(:,:)=(Energy(:,:,1)+Energy(:,:,2)+Energy(:,:,3)+Energy(:,:,4))/4;
homogeneity(:,:)=(Homogeneity(:,:,1)+Homogeneity(:,:,2)+Homogeneity(:,:,3)+Homogeneity(:,:,4))/4;
entropy(:,:)=(Entropy(:,:,1)+Entropy(:,:,2)+Entropy(:,:,3)+Entropy(:,:,4))/4;
end


        

