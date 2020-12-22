%%%%%%%%%%% The code for the EDI features extraction (before the clustering)
%%%%% Writed by Lan He
%%% Date: 2019-12-06
%%% E-mail: helan0811@126.com
%%% Department of Radiology, Guangdong Provincial People¡¯s Hospital, Guangdong Academy of Medical Sciences, Guangzhou 510080, China. 


clear all;
root=uigetdir;                %% read the file path for the file folder which contian the CT images and the VOI files                                                                          

addpath('NIfTI_20140122');    %% Adding the file folder of "NIFI_20140122"      
roiroot=cellstr([char(root),'/','roi.nii.gz']);            %% VOI files
imgroot=cellstr([char(root),'/','img.nii.gz']);            %% CT image files
roi=load_nii(char(roiroot));
ROI=roi.img;
[n1,n2,n3]=size(ROI);
A=[];
n=1;
for i=1:n3
    if sum(sum(ROI(:,:,i)))~=0
        A(n)=i;
        n=n+1;
    end
end


img=load_nii(char(imgroot));
IMG=img.img;
IMG=double(IMG);
ROI=double(ROI);
sd_total=zeros(1);
Mean_total=zeros(1);
Skewness_total=zeros(1);
Kurtosis_total=zeros(1);
contrast_total=zeros(1);
correlation_total=zeros(1);
energy_total=zeros(1);
homogeneity_total=zeros(1);
entropy_total=zeros(1);
num=1;


%%%% no filter
for i=1:length(A)
    im=IMG(:,:,A(i));
    roimarsk=ROI(:,:,A(i));
    msize=size(im);
    [I1,I2]=find(ROI(:,:,A(i)));
    roilength=length(I1);
    im1=zeros(msize(1),msize(2));
    im1=im.*roimarsk;
    im2=imcrop(im1,[min(I2)-2 min(I1)-2 max(I2)-min(I2)+3 max(I1)-min(I1)+3]);
    [sd,Mean,Skewness, Kurtosis]=extract_his(im2);
    [contrast,correlation,energy,homogeneity,entropy]=extract_GLCM(im2);
    [h,w]=size(sd);
    for i=1:h
        for j=1:w
            if ~(sd(i,j))==0
            sd_total(num)=sd(i,j);
            Mean_total(num)=Mean(i,j);
            Skewness_total(num)=Skewness(i,j);
            Kurtosis_total(num)=Kurtosis(i,j);
            contrast_total(num)=contrast(i,j);
            correlation_total(num)=correlation(i,j);
            energy_total(num)=energy(i,j);
            homogeneity_total(num)=homogeneity(i,j);
            entropy_total(num)=entropy(i,j);
            num=num+1;
            end
        end
    end
    
end
l=length(sd_total);
f=zeros(9,l);
f(1,:)=roundn(sd_total,-1);   
f(2,:)=roundn(Mean_total,-2);
f(3,:)=roundn(Skewness_total,-2);
f(4,:)=roundn(Kurtosis_total,-2);
f(5,:)=roundn(contrast_total,-2);
f(6,:)=roundn(correlation_total,-2);
f(7,:)=roundn(energy_total,-2);
f(8,:)=roundn(homogeneity_total,-2);
f(9,:)=roundn(entropy_total,-2);

name1=char(root);
name1=char(cellstr([name1,'_0.csv']));          %%write the result into the excel table
csvwrite(name1,f);

%%%% filter=1.0
sd_total=zeros(1);
Mean_total=zeros(1);
Skewness_total=zeros(1);
Kurtosis_total=zeros(1);
contrast_total=zeros(1);
correlation_total=zeros(1);
energy_total=zeros(1);
homogeneity_total=zeros(1);
entropy_total=zeros(1);
num=1;
for i=1:length(A)
    im=IMG(:,:,A(i));
    roimarsk=ROI(:,:,A(i));
    imf=lapFilter(im,5,1.0);
    msize=size(im);
    [I1,I2]=find(ROI(:,:,A(i)));
    roilength=length(I1);
    im1=zeros(msize(1),msize(2));
    im1=imf.*roimarsk;
    im2=imcrop(im1,[min(I2)-2 min(I1)-2 max(I2)-min(I2)+3 max(I1)-min(I1)+3]);
    [sd,Mean,Skewness, Kurtosis]=extract_his(im2);
    [contrast,correlation,energy,homogeneity,entropy]=extract_GLCM(im2);
    [h,w]=size(sd);
    for i=1:h
        for j=1:w
            if ~(sd(i,j))==0
            sd_total(num)=sd(i,j);
            Mean_total(num)=Mean(i,j);
            Skewness_total(num)=Skewness(i,j);
            Kurtosis_total(num)=Kurtosis(i,j);
            contrast_total(num)=contrast(i,j);
            correlation_total(num)=correlation(i,j);
            energy_total(num)=energy(i,j);
            homogeneity_total(num)=homogeneity(i,j);
            entropy_total(num)=entropy(i,j);
            num=num+1;
            end
        end
    end
    
end
l=length(sd_total);
f=zeros(9,l);
f(1,:)=roundn(sd_total,-1);   
f(2,:)=roundn(Mean_total,-2);
f(3,:)=roundn(Skewness_total,-2);
f(4,:)=roundn(Kurtosis_total,-2);
f(5,:)=roundn(contrast_total,-2);
f(6,:)=roundn(correlation_total,-2);
f(7,:)=roundn(energy_total,-2);
f(8,:)=roundn(homogeneity_total,-2);
f(9,:)=roundn(entropy_total,-2);

name1=char(root);
name1=char(cellstr([name1,'_1.0.csv']));          %%write the result into the excel table
csvwrite(name1,f);


%%%% filter=1.5
sd_total=zeros(1);
Mean_total=zeros(1);
Skewness_total=zeros(1);
Kurtosis_total=zeros(1);
contrast_total=zeros(1);
correlation_total=zeros(1);
energy_total=zeros(1);
homogeneity_total=zeros(1);
entropy_total=zeros(1);
num=1;
for i=1:length(A)
    im=IMG(:,:,A(i));
    roimarsk=ROI(:,:,A(i));
    imf=lapFilter(im,5,1.5);
    msize=size(im);
    [I1,I2]=find(ROI(:,:,A(i)));
    roilength=length(I1);
    im1=zeros(msize(1),msize(2));
    im1=imf.*roimarsk;
    im2=imcrop(im1,[min(I2)-2 min(I1)-2 max(I2)-min(I2)+3 max(I1)-min(I1)+3]);
    [sd,Mean,Skewness, Kurtosis]=extract_his(im2);
    [contrast,correlation,energy,homogeneity,entropy]=extract_GLCM(im2);
    [h,w]=size(sd);
    for i=1:h
        for j=1:w
            if ~(sd(i,j))==0
            sd_total(num)=sd(i,j);
            Mean_total(num)=Mean(i,j);
            Skewness_total(num)=Skewness(i,j);
            Kurtosis_total(num)=Kurtosis(i,j);
            contrast_total(num)=contrast(i,j);
            correlation_total(num)=correlation(i,j);
            energy_total(num)=energy(i,j);
            homogeneity_total(num)=homogeneity(i,j);
            entropy_total(num)=entropy(i,j);
            num=num+1;
            end
        end
    end
    
end

l=length(sd_total);
f=zeros(9,l);
f(1,:)=roundn(sd_total,-1);   
f(2,:)=roundn(Mean_total,-2);
f(3,:)=roundn(Skewness_total,-2);
f(4,:)=roundn(Kurtosis_total,-2);
f(5,:)=roundn(contrast_total,-2);
f(6,:)=roundn(correlation_total,-2);
f(7,:)=roundn(energy_total,-2);
f(8,:)=roundn(homogeneity_total,-2);
f(9,:)=roundn(entropy_total,-2);

name1=char(root);
name1=char(cellstr([name1,'_1.5.csv']));          %%write the result into the excel table
csvwrite(name1,f);


%%%% filter=2.0
sd_total=zeros(1);
Mean_total=zeros(1);
Skewness_total=zeros(1);
Kurtosis_total=zeros(1);
contrast_total=zeros(1);
correlation_total=zeros(1);
energy_total=zeros(1);
homogeneity_total=zeros(1);
entropy_total=zeros(1);
num=1;
for i=1:length(A)
    im=IMG(:,:,A(i));
    roimarsk=ROI(:,:,A(i));
    imf=lapFilter(im,5,2.0);
    msize=size(im);
    [I1,I2]=find(ROI(:,:,A(i)));
    roilength=length(I1);
    im1=zeros(msize(1),msize(2));
    im1=imf.*roimarsk;
    im2=imcrop(im1,[min(I2)-2 min(I1)-2 max(I2)-min(I2)+3 max(I1)-min(I1)+3]);
    [sd,Mean,Skewness, Kurtosis]=extract_his(im2);
    [contrast,correlation,energy,homogeneity,entropy]=extract_GLCM(im2);
    [h,w]=size(sd);
    for i=1:h
        for j=1:w
            if ~(sd(i,j))==0
            sd_total(num)=sd(i,j);
            Mean_total(num)=Mean(i,j);
            Skewness_total(num)=Skewness(i,j);
            Kurtosis_total(num)=Kurtosis(i,j);
            contrast_total(num)=contrast(i,j);
            correlation_total(num)=correlation(i,j);
            energy_total(num)=energy(i,j);
            homogeneity_total(num)=homogeneity(i,j);
            entropy_total(num)=entropy(i,j);
            num=num+1;
            end
        end
    end
    
end

l=length(sd_total);
f=zeros(9,l);
f(1,:)=roundn(sd_total,-1);   
f(2,:)=roundn(Mean_total,-2);
f(3,:)=roundn(Skewness_total,-2);
f(4,:)=roundn(Kurtosis_total,-2);
f(5,:)=roundn(contrast_total,-2);
f(6,:)=roundn(correlation_total,-2);
f(7,:)=roundn(energy_total,-2);
f(8,:)=roundn(homogeneity_total,-2);
f(9,:)=roundn(entropy_total,-2);

name1=char(root);
name1=char(cellstr([name1,'_2.0.csv']));          %%write the result into the excel table
csvwrite(name1,f);


%%%% filter=2.5
sd_total=zeros(1);
Mean_total=zeros(1);
Skewness_total=zeros(1);
Kurtosis_total=zeros(1);
contrast_total=zeros(1);
correlation_total=zeros(1);
energy_total=zeros(1);
homogeneity_total=zeros(1);
entropy_total=zeros(1);
num=1;
for i=1:length(A)
    im=IMG(:,:,A(i));
    roimarsk=ROI(:,:,A(i));
    imf=lapFilter(im,5,2.5);
    msize=size(im);
    [I1,I2]=find(ROI(:,:,A(i)));
    roilength=length(I1);
    im1=zeros(msize(1),msize(2));
    im1=imf.*roimarsk;
    im2=imcrop(im1,[min(I2)-2 min(I1)-2 max(I2)-min(I2)+3 max(I1)-min(I1)+3]);
    [sd,Mean,Skewness, Kurtosis]=extract_his(im2);
    [contrast,correlation,energy,homogeneity,entropy]=extract_GLCM(im2);
    [h,w]=size(sd);
    for i=1:h
        for j=1:w
            if ~(sd(i,j))==0
            sd_total(num)=sd(i,j);
            Mean_total(num)=Mean(i,j);
            Skewness_total(num)=Skewness(i,j);
            Kurtosis_total(num)=Kurtosis(i,j);
            contrast_total(num)=contrast(i,j);
            correlation_total(num)=correlation(i,j);
            energy_total(num)=energy(i,j);
            homogeneity_total(num)=homogeneity(i,j);
            entropy_total(num)=entropy(i,j);
            num=num+1;
            end
        end
    end
    
end

l=length(sd_total);
f=zeros(9,l);
f(1,:)=roundn(sd_total,-1);   
f(2,:)=roundn(Mean_total,-2);
f(3,:)=roundn(Skewness_total,-2);
f(4,:)=roundn(Kurtosis_total,-2);
f(5,:)=roundn(contrast_total,-2);
f(6,:)=roundn(correlation_total,-2);
f(7,:)=roundn(energy_total,-2);
f(8,:)=roundn(homogeneity_total,-2);
f(9,:)=roundn(entropy_total,-2);

name1=char(root);
name1=char(cellstr([name1,'_2.5.csv']));          %%write the result into the excel table
csvwrite(name1,f);
