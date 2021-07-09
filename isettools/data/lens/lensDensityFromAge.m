function TL = lensDensityFromAge(A,wave)
% Compute the spectral density of the lens at different ages
%
%  Synopsis
%    TL = lensDensityFromAge(AgeInYears,wave)
%
% See also
%   lensPokornySmith.m

% Example:
%{
wave = 400:1:700;
ages = [10:5:70];
lensDensity = lensDensityFromAge(ages,wave);
mesh(ages,wave,lensDensity);
ieNewGraphWin; plot(wave,lensDensity); grid on; xlabel('Wavelength');
%}

%%
if ~exist('wave','var'), wave = 400:10:700; end

d = ieReadSpectra('lensPokornySmith',wave);
TL1 = d(:,1);
TL2 = d(:,2);

%% Formula for Pokorny and Smith, but use lensDensity default for younger ages

nAges = numel(A);
TL = zeros(numel(wave),nAges);

for ii=1:nAges
    if A(ii) <= 20
        TL(:,ii) = ieReadSpectra('lensDensity.mat',wave);
    elseif 20 < A(ii) && A(ii) < 60
        TL(:,ii) = TL1*(1 + 0.02 * (A(ii) - 32)) + TL2;
        
    elseif A(ii) >= 60
        TL(:,ii) = TL1*(1.56 + 0.0667*(A(ii) - 60)) + TL2;
    end
    
end

end

%%

