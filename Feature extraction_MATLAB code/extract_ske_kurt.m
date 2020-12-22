function [f_skewness_mean,f_kurtosis_mean]=extract_ske_kurt(im,roimask,patch_size)

if mod(patch_size,2) == 0                           %mod为取模，判断是否是2的倍数                                 
    patch_size = patch_size + 1;
end
border_size = (patch_size - 1)/2;
[height,width]=size(im);
f_skewness = zeros(height, width);
f_kurtosis = zeros(height, width);
blocks = im2col(im,[patch_size patch_size],'sliding');
%--------------------------------------- skewness -------------------------------------------------------%     
    f = col2im(skewness(blocks),[patch_size patch_size],[height, width],'sliding');
    f(isnan(f)) = 0;
    f_skewness(border_size+1:end-border_size, border_size+1:end-border_size) = f; 
%--------------------------------------- kurtosis -------------------------------------------------------%    
    f = col2im(kurtosis(blocks),[patch_size patch_size],[height, width],'sliding');
    f(isnan(f)) = 0;
    f_kurtosis(border_size+1:end-border_size, border_size+1:end-border_size) = f;
%---------------------------------------------------------------------------------------------------------%
f_skewness_mean=mean(f_skewness(roimask));
f_kurtosis_mean=mean(f_kurtosis(roimask));
end

