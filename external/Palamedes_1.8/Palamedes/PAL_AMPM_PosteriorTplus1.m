%
%PAL_AMPM_PosteriorTplus1  Derives posterior distributions for combinations
%   of stimulus intensity and response on trial T + 1 in psi method 
%   adaptive procedure.
%   
%   syntax: [PosteriorTplus1givenSuccess PosteriorTplus1givenFailure ...
%       pSuccessGivenx] = PAL_AMPM_PosteriorTplus1(pdf, PFLookUpTable)
%
%   Internal function
%
%Introduced: Palamedes version 1.0.0 (NP)
%Modified: Palamedes version 1.5.0, 1.6.1, 1.6.3 (see History.m)

function [PosteriorTplus1givenSuccess, PosteriorTplus1givenFailure, pSuccessGivenx] = PAL_AMPM_PosteriorTplus1(pdf, PFLookUpTable)

pdf5D = repmat(pdf, [1 1 1 1 size(PFLookUpTable,5)]);
pSuccessGivenx = sum(sum(sum(sum(pdf5D.*PFLookUpTable,1),2),3),4);

PosteriorTplus1givenSuccess = pdf5D.*PFLookUpTable;
PosteriorTplus1givenFailure = pdf5D-PosteriorTplus1givenSuccess;

if exist('bsxfun.m','file')
    PosteriorTplus1givenSuccess = bsxfun(@rdivide,PosteriorTplus1givenSuccess,sum(sum(sum(sum(PosteriorTplus1givenSuccess,1),2),3),4));
    PosteriorTplus1givenFailure = bsxfun(@rdivide,PosteriorTplus1givenFailure,sum(sum(sum(sum(PosteriorTplus1givenFailure,1),2),3),4));
else    
    Denominator = squeeze(sum(sum(sum(sum(PosteriorTplus1givenSuccess,1),2),3),4));
    Denominator = repmat(Denominator, [1 size(pdf5D,1) size(pdf5D,2) size(pdf5D,3) size(pdf5D,4)]);
    Denominator = permute(Denominator, [2 3 4 5 1]);
    PosteriorTplus1givenSuccess = PosteriorTplus1givenSuccess./Denominator;
    
    Denominator = squeeze(sum(sum(sum(sum(PosteriorTplus1givenFailure,1),2),3),4));
    Denominator = repmat(Denominator, [1 size(pdf5D,1) size(pdf5D,2) size(pdf5D,3) size(pdf5D,4)]);
    Denominator = permute(Denominator, [2 3 4 5 1]);
    PosteriorTplus1givenFailure = PosteriorTplus1givenFailure./Denominator;
end