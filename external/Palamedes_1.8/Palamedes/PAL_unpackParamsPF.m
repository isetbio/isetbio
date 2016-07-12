%
%PAL_unpackParamsPF     Re-format parameters of PF input
%
%syntax: [alpha beta gamma lambda] = PAL_unpackParamsPF(params)
%
%Internal function
%
%Introduced: Palamedes version 1.0.0 (NP)
%Modified: Palamedes version 1.6.3 (see History.m)

function [alpha, beta, gamma, lambda] = PAL_unpackParamsPF(params)

gamma = 0;
lambda = 0;

if isstruct(params)
    alpha = params.alpha;
    beta = params.beta;
    if isfield(params,'gamma')
        gamma = params.gamma;
        if isfield(params,'lambda')
            lambda = params.lambda;
        end
    end
else
    alpha = params(1);
    beta = params(2);
    if length(params) > 2
        gamma = params(3);
        if length(params) > 3
            lambda = params(4);
        end
    end
end