%-----------------------------------------------------------------------------
function E = calculateEntropy(glcm)
glcm(glcm == 0) = [];
foo = (-1)*glcm.*log(glcm);
E = sum(foo(:));
end