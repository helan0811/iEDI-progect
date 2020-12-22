function im_filt = lapFilter(im, win, sigma)
% Zhang Zhongping
if mod(win,2) == 0
    win = win+1;
end
hwin = (win-1) / 2;

% Laplacian Gaussian filter
filter = zeros(win);
for i = -hwin : hwin
    for j = -hwin : hwin
        exponent = -(i.^2+j.^2) / (2*sigma.^2);
        filter(i+hwin+1, j+hwin+1) = 1 / (pi*sigma.^4) * (1+exponent)*exp(exponent);
        % whether should there be -1
        %filter(i+hwin+1, j+hwin+1) = -1 / (pi*sigma.^4) * (1+exponent)*exp(exponent);
    end
end

im_filt = im;
[nrow ncol] = size(im);

for i = 1+hwin:ncol-hwin
    for j = 1+hwin:nrow-hwin               
        tmp = im(i-hwin:i+hwin, j-hwin:j+hwin) .* filter;
        im_filt(i, j) = sum(tmp(:));     
    end
end

end