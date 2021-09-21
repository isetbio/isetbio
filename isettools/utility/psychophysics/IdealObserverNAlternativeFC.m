function pCorrect = IdealObserverNAlternativeFC(meanResponses,nSimulatedTrials)
% Monte-Carlo simulation of probability correct for N-alternative forced choice
%
% Synopsis:
%    probCorrect = IdealObserverNAlternativeFC(meanResponses,nSimulatedTrials)
%
% Description:
%
% Inputs:
%    meanReponses -        Matrix (nResponses x nAlternatives).  Each
%                          is the vector of mean sensor responses (e.g.
%                          excitations of the cones in the mosaic) for each
%                          of the alternatives.  Actual responses are
%                          assumed to be Poisson distributed with these
%                          means.
%    nSimulatedTrials -    Number of Monte-Carlo simulated trials to use.
%
% Outputs:
%    probCorrect -         Probability correct
%
% Optional Key-Value Pairs

% History
%    09/21/21   dhb  Wrote it.

% Examples:
%{
    % Computes probability correct versus delta for 4AFC, where the alternatives
    % have mean vectors with all entries but 1 equal to base, and with one entry
    % equal to base plus delta.  Each alternative has a different entry
    % perturbed by delta.  As expected for Poisson noise, the task gets
    % harder for fixed delta as base increases, and easier for fixed base
    % delta increasess.  The plot shows the latter effect.  Vary the value
    % of base to see the former.  Using nSimulatedTrials = 1000 gives a
    % quick answer that is a little noisy - increase for more precision.
    clear; close all;
    base = 10;
    deltas = [0 1 2 5 10 20 40 80];
    nSimulatedTrials = 1000;
    for i = 1:length(deltas)
    meanResponses = [ [base+deltas(i) base base base]', [base base+deltas(i) base base]', [base base base+deltas(i) base]', [base base base base+deltas(i) ]'];
    pCorrects(i) = IdealObserverNAlternativeFC(meanResponses,nSimulatedTrials);
    fprintf('Four alternative FC, base = %d, delta = %d, pCorrect = %0.2g\n', base, deltas(i), pCorrects(i));
    end
    figure; clf; hold on;
    plot(deltas,pCorrects,'ro','MarkerFaceColor','r','MarkerSize',12);
    plot(deltas,pCorrects,'r');
    xlabel('Delta');
    ylabel('pCorrect');
    xlim([0 100]); ylim([0 1]);
%}


% Initialize
nCorrect = 0;
[nResponses,nAlternatives] = size(meanResponses);

% Simulated trials
for ii = 1:nSimulatedTrials
    % Chose alternative
    whichAlternative = randi(nAlternatives);

    % Get responses
    theResponses = poissrnd(meanResponses(:,whichAlternative));

    % Compute log likelihood of each response
    logLikelihoods = zeros(1,nAlternatives);
    for jj = 1:nAlternatives
        logLikelihoods(jj) = PoissonLogLikelihoodCompare(theResponses,meanResponses(:,jj));
    end

    % Find max
    [~,maxIndex] = max(logLikelihoods);

    % Was it right?
    if (maxIndex == whichAlternative)
        nCorrect = nCorrect + 1;
    end
end

% Compute probability correct
pCorrect = nCorrect/nSimulatedTrials;

end

% Compute the part of the Poisson log likelihood that depends on the mean
% responses.  It's very slow to the term ln(theRespones!), which is needed for
% the full log likelihood, but we only need to compare across the various
% choices of mean responses, and thus don't need this term.
%
% To derive this for each entry, take the log of the probability mass function 
% and drop the term -ln(theResponses!).  The note that the joint liklihood of
% all the entries is the product of the individual entry likelihoods and that
% the log turns this into the sum of the individual log likihoods.
function logLikelihood = PoissonLogLikelihoodCompare(theResponses,meanResponses)

logLikelihood = sum(-meanResponses + log(meanResponses).*theResponses);

end