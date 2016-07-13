%
%PAL_PFBA_Fit   Fit threshold and slope parameters of psychometric function 
%   using Bayesian criterion.
%
%syntax: [paramsValues posterior] = PAL_PFBA_Fit(StimLevels, NumPos, ...
%           OutOfNum, priorAlphaValues, priorBetaValues, gamma, lambda, ...
%           PF, {optional arguments})
%
%PAL_PFBA_Fit derives a posterior distribution across threshold x slope
%   parameter space and computes expected value of threshold and slope (to
%   be used as parameter estimates). Standard errors are derived as the
%   marginal standard deviations of threshold and slope in posterior
%   distribution.
%
%Input: 
%   'StimLevels': vector containing the stimulus levels utilized
%
%   'NumPos': vector of equal length to 'StimLevels' containing for each of 
%       the stimulus levels the number of positive responses (e.g., ‘yes’ or 
%       ‘correct’) observed.
%
%   'OutOfNum': vector of equal length to 'StimLevels' containing the 
%       number of trials tested at each of the stimulus levels.
%
%   'priorAlphaValues': vector of any length which specifies which 
%       threshold values are to be included in the prior's parameter space.
%
%   'priorBetaValues': vector of any length which specifies which slope 
%       values are to be included in the prior's parameter space.
%
%   'gamma': scalar corresponding to the guess rate to be assumed.
%
%   'lambda': scalar corresponding to the lapse rate to be assumed.
%
%   'PF': The psychometric function to be fitted. This needs to be passed 
%       as an inline function. Options include:
%           @PAL_Logistic
%           @PAL_Weibull
%           @PAL_Gumbel (i.e., log-Weibull)
%           @PAL_Quick
%           @PAL_logQuick
%           @PAL_CumulativeNormal
%           @PAL_Gumbel
%           @PAL_HyperbolicSecant
%
%Output:
%   'paramsValues': four-element vector containing threshold estimate,
%       slope estimate, standard error of threshold estimate, and standard
%       error of slope estimate respectively.
%
%   'posterior': The posterior distribution. Useful to check e.g., whether 
%       the posterior is (effectively) contained within the parameter space
%       considered or whether the distribution is symmetrical.
%
%   User has the option to specify a prior distribution (if no prior is
%       provided, a rectangular distribution on parameter space is used as 
%       prior [aka 'uniform prior' or 'prior of ignorance']). The prior 
%       distribution should be defined across the parameter space [i.e., be 
%       of size length(priorBetaValues)xlength(priorAlphaValues)]. Use 
%       optional argument 'prior', followed by array containing prior to 
%       specify a prior other than the uniform prior. 
%
%Full Example:
%
%   StimLevels = [-3:1:3];
%   NumPos = [55 55 66 75 91 94 97];
%   OutOfNum = [100 100 100 100 100 100 100];
%   priorAlphaValues = [-1:.01:1];
%   priorBetaValues = [-.5:.01:.5];
%   gamma = 0.5;
%   lambda = 0;
%   PF = @PAL_Logistic;
%
%   [Alpha Beta] = meshgrid(priorAlphaValues,priorBetaValues);
%   prior = PAL_pdfNormal(Alpha,0,1).*PAL_pdfNormal(Beta,0,1);
%   prior = prior./sum(sum(prior)); %such that sum(sum(Prior)) == 1
%
%   [paramsValues posterior] = PAL_PFBA_Fit(StimLevels, NumPos, ...
%       OutOfNum, priorAlphaValues, priorBetaValues, gamma, lambda, PF, ...
%       'prior', prior);
%
%Introduced: Palamedes version 1.0.0 (NP)
%Modified: Palamedes version 1.6.3 (see History.m)

function [paramsValues, posterior] = PAL_PFBA_Fit(StimLevels, NumPos, OutOfNum, priorAlphaValues, priorBetaValues, gamma, lambda, PF, varargin)

prior = ones(length(priorBetaValues),length(priorAlphaValues))/(length(priorBetaValues).*length(priorAlphaValues));

if ~isempty(varargin)
    NumOpts = length(varargin);
    for n = 1:2:NumOpts
        valid = 0;
        if strcmpi(varargin{n}, 'prior')
            prior = varargin{n+1};
            valid = 1;
        end
        if valid == 0
            warning('PALAMEDES:invalidOption','%s is not a valid option. Ignored.',varargin{n});
        end        
    end            
end

prior = log(prior);

[StimLevels, NumPos, OutOfNum] = PAL_PFML_GroupTrialsbyX(StimLevels, NumPos, OutOfNum);

[params.alpha, params.beta] = meshgrid(priorAlphaValues,priorBetaValues);
params.beta = 10.^params.beta;
params.gamma = gamma;
params.lambda = lambda;

llikelihood = prior;

for StimLevel = 1:length(StimLevels)
    pcorrect = PF(params, StimLevels(StimLevel));
    llikelihood = llikelihood+log(pcorrect).*NumPos(StimLevel)+log(1-pcorrect).*(OutOfNum(StimLevel)-NumPos(StimLevel));
end

llikelihood = llikelihood-max(max(llikelihood)); %this avoids likelihoods containing zeros
likelihood = exp(llikelihood);
posterior = likelihood./sum(sum(likelihood));

posteriora = sum(posterior,1);
posteriorb = sum(posterior,2);

paramsValues(1) = sum(posteriora.*priorAlphaValues);
paramsValues(2) = sum(posteriorb.*priorBetaValues');
paramsValues(3) = (sum(posteriora.*(priorAlphaValues-paramsValues(1)).^2))^.5;
paramsValues(4) = (sum(posteriorb.*(priorBetaValues'-paramsValues(2)).^2))^.5;