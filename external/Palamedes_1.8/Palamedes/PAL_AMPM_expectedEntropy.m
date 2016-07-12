%
%PAL_AMPM_expectedEntropy  expected entropy values for different stimulus
%   intensities.
%
%   syntax: [PM expectedEntropy] = PAL_AMPM_expectedEntropy(PM,varargin)
%
%Internal function
%
%Introduced: Palamedes version 1.6.0 (NP)
%Modified: Palamedes version 1.6.1, 1.6.3 (see History.m)

function [PM, expectedEntropy] = PAL_AMPM_expectedEntropy(PM,varargin)

fixLapse = false;

valid = 0;
if ~isempty(varargin)
    if strncmpi(varargin{1}, 'fixL',4)            
        fixLapse = varargin{2};                
        valid = 1;
    end
    if valid == 0
        warning('PALAMEDES:invalidOption','%s is not a valid option. Ignored.',varargin{1})
    end        
end

[PM.posteriorTplus1givenSuccess, PM.posteriorTplus1givenFailure, pSuccessGivenx] = PAL_AMPM_PosteriorTplus1(PM.pdf, PM.LUT);

if fixLapse
    tempSuccess = zeros(size(PM.posteriorTplus1givenSuccess(:,:,:,1,:)));
    tempFailure = zeros(size(PM.posteriorTplus1givenSuccess(:,:,:,1,:)));
    [trash, I] = PAL_findMax(PM.pdf);
    for SR = 1:length(PM.stimRange)
        tempSuccess(:,:,:,:,SR) = PM.posteriorTplus1givenSuccess(:,:,:,I(4),SR)./sum(sum(sum(sum(sum(PM.posteriorTplus1givenSuccess(:,:,:,I(4),SR))))));
        tempFailure(:,:,:,:,SR) = PM.posteriorTplus1givenFailure(:,:,:,I(4),SR)./sum(sum(sum(sum(sum(PM.posteriorTplus1givenFailure(:,:,:,I(4),SR))))));
    end
    expectedEntropy = PAL_Entropy(tempSuccess,4,PM.marginalize).*pSuccessGivenx + PAL_Entropy(tempFailure,4,PM.marginalize).*(1-pSuccessGivenx);        
else
    expectedEntropy = PAL_Entropy(PM.posteriorTplus1givenSuccess,4,PM.marginalize).*pSuccessGivenx + PAL_Entropy(PM.posteriorTplus1givenFailure,4,PM.marginalize).*(1-pSuccessGivenx);        
end