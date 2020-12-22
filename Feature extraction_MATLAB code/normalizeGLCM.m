%-----------------------------------------------------------------------------
function glcm = normalizeGLCM(glcm)
% Normalize glcm so that sum(glcm(:)) is one.
if any(glcm(:))
  glcm = glcm ./ sum(glcm(:));
end
end