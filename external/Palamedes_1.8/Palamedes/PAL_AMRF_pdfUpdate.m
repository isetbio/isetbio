%
%PAL_AMRF_pdfUpdate Updates posterior distribution in running fit adaptive
%   method based on response to a trial.
%
%   syntax: [pdf] = PAL_AMRF_pdfUpdate(pdf, alphas, beta, gamma, ...
%       lambda, x, response, PF)
%   
%   Internal function
%
%Introduced: Palamedes version 1.0.0 (NP)


function [pdf] = PAL_AMRF_pdfUpdate(pdf, alphas, beta, gamma, lambda, x, response, PF)

params.alpha = alphas;
params.beta = beta;
params.gamma = gamma;
params.lambda = lambda;

p = PF(params,x);

if response == 1
    pdf = pdf.*p;
else
    pdf = pdf.*(1 - p);
end

pdf = pdf./sum(pdf);