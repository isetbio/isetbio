function lsf = psf2lsf(psf)
% Derive the line spread from the pointspread
%
% Input
%   psf - Pointspread function.  3rd dimension is wavelength
%
% Optional key/val
%    N/A
%
% Output
%    lsf - line spread function.  The columns are the lsf for each wave
%
% See also
%    psf2otf

% Examples:
%{
   psf = randn(128,128);
   lsf = psf2lsf(psf);
%}
%{
   psf = randn(10,10,3);
   lsf = psf2lsf(psf);
%}

nWave = size(psf,3);

if nWave > 1
    lsf = zeros(size(psf,2),nWave);
    for ii = 1:nWave
        lsf(:,ii) = sum(psf(:,:,ii),1);
    end
else
    lsf = sum(psf,1);
end

end
