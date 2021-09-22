function [pCorrect,predCorrect,whichAlternatives,decisionVariables, ...
    meanDecisionGiven1,meanDecisionGiven2, ...
    predMeanDecisionGiven1,predMeanDecisionGiven2,predVarDecisionGiven1,predVarDecisionGiven2, ...
    logLikelihoodsCompare1,logLikelihoodsCompare2] = PoissonIdealObserverNAlternativeFC(meanResponses,nSimulatedTrials,DEBUG)
% Monte-Carlo simulation of probability correct for N-alternative forced choice
%
% Synopsis:
%    [pCorrect,whichAlternatives,decisionVariables, ...
%    meanDecisionGiven1,meanDecisionGiven2, ...
%    predMeanDecisionGiven1,predMeanDecisionGiven2,predVarDecisionGiven1,predVarDecisionGiven2] = PoissonIdealObserverNAlternativeFC(meanResponses,nSimulatedTrials)
%
% Description:
%    Simulated performance on an N alternative forced choice task.  See
%    examples in source code for how to call to model a 4AFC task with one
%    stimulus presented on each trial and the subject indicating which one,
%    and a TAFC task where both stimuli are presented on each trial.
%
% Inputs:
%    meanReponses     -    Matrix (nResponses x nAlternatives).  Each
%                          is the vector of mean sensor responses (e.g.
%                          excitations of the cones in the mosaic) for each
%                          of the alternatives.  Actual responses are
%                          assumed to be Poisson distributed with these
%                          means.
%    nSimulatedTrials -    Number of Monte-Carlo simulated trials to use.
%    DEBUG            -    Boolean (default false).  Returns extra
%                          information to try to understand things better.
%
% Outputs:
%    probCorrect -         Probability correct. 
%
% Outputs if nAlternatives == 2 and DEBUG is true.  These variables are set
% and returned even if these conditions are not meant, but are not meaningful
% in that case. Use of these extra returns is illustrated in the TAFC
% example in the source code.
%    predCorrect       -   Prediction of pCorrect based on normal
%                          approximation to the decision variable
%                          histogram.
%    whichAlternatives -   nSimulatedTrials vector containing 1 or 2
%                          according to which stimulus was simulated on each trial.
%    decisionVariables -   nSimulatedTrials vector containing the difference between
%                          the altnerative 2 and alternative 1 log
%                          likelihoods. Choose alternative 2 if this is
%                          greater than 0 and alternative 1 if less.
%    meanDecisionGiven1 -  Mean of decisionVariables when stimulus 1
%                          presented.
%    meanDecisionGiven2 -  Mean of decisionVariables when stimulus 2
%                          presented.
%    predMeanDecisionGiven1 - Analytic prediction of stimulus 1 mean above.
%    predMeanDecisionGiven2 - Analytic prediction of stimulus 2 mean above.
%    predVarDecisionGiven1 - Analytic prediction of variance of decision variable
%                          given sttimulus 1 presented. 
%    predVarDecisionGiven2 - Analytic prediction of variance of decision variable
%                          given sttimulus 1 presented. 
%
% Optional Key-Value Pairs
%    None.
%
% See also: analyticPoissonIdealObserver, SupportVectorMachineObserverNAlternativeFC

% History
%    09/21/21   dhb  Wrote it.

% Examples:
%{
    % Computes probability correct versus delta for one stimulus per trial
    % 4AFC, where the alternatives have mean vectors with all entries but 1
    % equal to base, and with one entry equal to base plus delta.  Each
    % alternative has a different entry perturbed by delta.  As expected
    % for Poisson noise, the task gets harder for fixed delta as base
    % increases, and easier for fixed base delta increasess.  The plot
    % shows the latter effect.  Vary the value of base to see the former.
    % Using nSimulatedTrials = 1000 gives a quick answer that is a little
    % noisy - increase for more precision.
    clear; close all;
    base = 100;
    deltas = [0 1 2 5 10 20 40 80];
    nSimulatedTrials = 1000;
    for i = 1:length(deltas)
        meanResponses = [ [base+deltas(i) base base base]', [base base+deltas(i) base base]', [base base base+deltas(i) base]', [base base base base+deltas(i) ]'];
        pCorrects(i) = PoissonIdealObserverNAlternativeFC(meanResponses,nSimulatedTrials);
        fprintf('Four alternative FC, base = %d, delta = %d, pCorrect = %0.2g\n', base, deltas(i), pCorrects(i));
    end
    figure; clf; hold on;
    plot(deltas,pCorrects,'ro-','MarkerFaceColor','r','MarkerSize',12);
    xlabel('Delta');
    ylabel('pCorrect');
    xlim([0 100]); ylim([0 1]);
    title(sprintf('4AFC Poisson Ideal Observer, base = %d',base));
%}
%{
    % Computes probability correct versus delta for 2AFC, with constant
    % addition differentiating the two means. This simulates two stimuli
    % per trial (e.g., blank and test in random order) and can be used
    % to compare to our analytic approximation for the same task.  The
    % way this is done is to pass the concatenated mean responses for the 
    % two stimuli, once in each order.
    %
    % This example compares to an analytic approximation from Geisler
    % (1984), see routine analyticPoissonIdealObserver.  Agreement is good
    % if you make nSimulatedTrials ~10000 or more. Also in agreement is 
    % prediction based on mean and variance of distribution of the decision
    % variable given each stimulus.
    %
    % This example also sets DEBUG to true and returns a bunch of extra
    % information that are useful if you want to delve into the details of
    % the analytic approximation in Geisler (1984). The example makes
    % various plots of this information.
    clear; close all;
    dimension = 10;
    base = 10;
    deltas = 0.04*[0 1 2 5 10 25 50 100];
    nSimulatedTrials = 1000;
    for i = 1:length(deltas)
        % Set up mean resopnses for each stimulus
        meanResponse1 = base*ones(dimension, 1);
        meanResponse2 = meanResponse1 + deltas(i)*ones(dimension,1);

        % Concatenate for TAFC where both stimuli are presented on each trial.
        meanResponses = [[meanResponse1 ; meanResponse2] , [meanResponse2 ; meanResponse1]];

        % If instead you wanted to model yes-no, set up mean responses as
        % commented out here. If you do this, the result will no longer
        % agree with what is computed in analyticPoissonIdealObserver.
        % meanResponses = [[meanResponse1] , [meanResponse2]];

        % Call this routine with DEBUG = true and collect up all of the
        % optional return values.
        [pCorrects(i), predCorrects(i), whichAlternatives{i},decisionVariables{i},meanDecisionsGiven1(i), meanDecisionsGiven2(i), ...
         predMeanDecisionsGiven1(i),predMeanDecisionsGiven2(i),predVarDecisionsGiven1(i),predVarDecisionsGiven2(i), ...
         logLikelihoodsCompare1{i}, logLikelihoodsCompare2{i}] = PoissonIdealObserverNAlternativeFC(meanResponses,nSimulatedTrials,true);

        % Call analyticaPoissonIdealObsever for comparison.
        [pCorrects1(i),nil,decisionMeans(i),decisionVars(i)] = analyticPoissonIdealObserver(meanResponse1,meanResponse2);
    end
    figure; clf; hold on;
    plot(deltas,pCorrects,'ro-','MarkerFaceColor','r','MarkerSize',12);
    plot(deltas,pCorrects1,'bo-','MarkerFaceColor','b','MarkerSize',10);
    plot(deltas,predCorrects,'go-','MarkerFaceColor','g','MarkerSize',8);
    xlabel('Delta');
    ylabel('pCorrect');
    xlim([0 max(deltas)]); ylim([0 1]);
    legend({'Monte Carlo','Approx (d-prime)','Approx (normcdf)'},'Location','SouthWest');
    title(sprintf('TAFC Poisson Ideal Observer, base = %d',base));

    % Histogram of decision variables given altnerative 1
    figure; clf; hold on;
    whichToPlot = 6;
    index = find(whichAlternatives{whichToPlot} == 1);
    h = histogram(decisionVariables{whichToPlot}(index),50,'Normalization','pdf');
    plot(h.BinEdges,normpdf(h.BinEdges,predMeanDecisionsGiven1(whichToPlot),sqrt(predVarDecisionsGiven1(whichToPlot))),'LineWidth',3);
    plot(mean(decisionVariables{whichToPlot}(index)),normpdf(mean(decisionVariables{whichToPlot}(index)),predMeanDecisionsGiven1(whichToPlot),sqrt(predVarDecisionsGiven1(whichToPlot))),'bo','MarkerFaceColor','b','MarkerSize',12); 

    % Histogram of decision variables given altnerative 2
    figure; clf; hold on;
    index = find(whichAlternatives{whichToPlot} == 2);
    h = histogram(decisionVariables{whichToPlot}(index),50,'Normalization','pdf');
    plot(h.BinEdges,normpdf(h.BinEdges,predMeanDecisionsGiven2(whichToPlot),sqrt(predVarDecisionsGiven2(whichToPlot))),'LineWidth',3);
    plot(mean(decisionVariables{whichToPlot}(index)),normpdf(mean(decisionVariables{whichToPlot}(index)),predMeanDecisionsGiven2(whichToPlot),sqrt(predVarDecisionsGiven2(whichToPlot))),'bo','MarkerFaceColor','b','MarkerSize',12); 

    % Plot of empirical and predicted decision variable means
    figure; clf; hold on;
    plot(predMeanDecisionsGiven1,meanDecisionsGiven1,'ro','MarkerFaceColor','r','MarkerSize',12);
    plot(predMeanDecisionsGiven2,meanDecisionsGiven2,'bo','MarkerFaceColor','b','MarkerSize',12);
    xlabel('Predicted Mean');
    ylabel('Monte Carlo Mean');
    limLow = min([predMeanDecisionsGiven1,meanDecisionsGiven1,predMeanDecisionsGiven2,meanDecisionsGiven2]);
    limHigh = max([predMeanDecisionsGiven1,meanDecisionsGiven1,predMeanDecisionsGiven2,meanDecisionsGiven2]);
    xlim([limLow limHigh]); ylim([limLow limHigh]);
    plot([limLow  limHigh],[limLow limHigh],'k');

    % Plot to see if log likelihoods under two hypotheses are correlated
    figure; clf; hold on;
    index = find(whichAlternatives{whichToPlot} == 1);
    plot(logLikelihoodsCompare2{whichToPlot}(index),logLikelihoodsCompare2{whichToPlot}(index)-logLikelihoodsCompare1{whichToPlot}(index),'ro','MarkerFaceColor','r','MarkerSize',12);
    index = find(whichAlternatives{whichToPlot} == 2);
    plot(logLikelihoodsCompare2{whichToPlot}(index),logLikelihoodsCompare2{whichToPlot}(index)-logLikelihoodsCompare1{whichToPlot}(index),'bo','MarkerFaceColor','b','MarkerSize',12);
%}

% Set default for DEBUG
if (nargin < 3 | isempty(DEBUG))
    DEBUG = false;
end

% Initialize
nCorrect = 0;
[nResponses,nAlternatives] = size(meanResponses);

% If nAlternatives is two, we return the decision variable
% as the difference between the two log likelihoods.  This
% is mainly here so that we can see how far the empirical
% distribution differs from a normal, as the normal approximation
% is used in an analytic approximation.
if (DEBUG & nAlternatives == 2)
    decisionVariables = zeros(nSimulatedTrials,1);
    sumDecisionGiven1 = 0;
    sumDecisionGiven2 = 0;
    n1 = 0;
    n2 = 0;
else
    decisionVariables = [];
end

% Simulated trials
for ii = 1:nSimulatedTrials
    % Chose alternative
    whichAlternative = randi(nAlternatives);

    % Get responses
    theResponses = poissrnd(meanResponses(:,whichAlternative));

    % Compute log likelihood of each response
    logLikelihoodCompare = zeros(1,nAlternatives);
    for jj = 1:nAlternatives
        logLikelihoodCompare(jj) = PoissonLogLikelihoodCompare(theResponses,meanResponses(:,jj));
    end

    % Find max
    [~,maxIndex] = max(logLikelihoodCompare);

    % Was it right?
    if (maxIndex == whichAlternative)
        nCorrect = nCorrect + 1;
    end

    % Accumulate decision variable if we are doing so.  Need DEBUG set to
    % true and this only is done if there are two alternatives.
    if (DEBUG & nAlternatives == 2)
        C = sum(meanResponses(:,1)-meanResponses(:,2));
        whichAlternatives(ii) = whichAlternative;
        decisionVariables(ii) = logLikelihoodCompare(2)-logLikelihoodCompare(1) ;
        logLikelihoodsCompare1(ii) = logLikelihoodCompare(1);
        logLikelihoodsCompare2(ii) = logLikelihoodCompare(2);
        if (whichAlternative == 1)
            sumDecisionGiven1 = sumDecisionGiven1 +  decisionVariables(ii);
            n1 = n1 + 1;
        else
            sumDecisionGiven2 = sumDecisionGiven2 +  decisionVariables(ii);
            n2 = n2 + 1;
        end
    end
end

% Compute probability correct
pCorrect = nCorrect/nSimulatedTrials;

% Mean of decision variables
if (DEBUG & nAlternatives == 2)
    meanDecisionGiven1 = sumDecisionGiven1/n1;
    meanDecisionGiven2 = sumDecisionGiven2/n2;
    predMeanDecisionGiven1 = sum(meanResponses(:,1) .* log(meanResponses(:,2)./meanResponses(:,1))) + C;
    predMeanDecisionGiven2 = sum(meanResponses(:,2) .* log(meanResponses(:,2)./meanResponses(:,1))) + C;
    predVarDecisionGiven1 = sum(meanResponses(:,1) .* (log(meanResponses(:,2)./meanResponses(:,1))).^2);
    predVarDecisionGiven2 = sum(meanResponses(:,2) .* (log(meanResponses(:,2)./meanResponses(:,1))).^2);
    % fprintf('Mean Z given 1: %g, predicted: %g\n',meanDecisionGiven1,predMeanDecisionGiven1);
    % fprintf('Mean Z given 2: %g, predicted: %g\n\n',meanDecisionGiven2,predMeanDecisionGiven2);

    % Given the normal approximation to the decision variable
    predCorrectGiven1 = normcdf(0,predMeanDecisionGiven1,sqrt(predVarDecisionGiven1));
    predCorrectGiven2 = 1-normcdf(0,predMeanDecisionGiven2,sqrt(predVarDecisionGiven2));
    predCorrect = (predCorrectGiven1 + predCorrectGiven2)/2;
else
    meanDecisionGiven1 = 0;
    meanDecisionGiven2 = 0;
    predMeanDecisionGiven1 = 0;
    predMeanDecisionGiven2 = 0;
    predVarDecisionGiven1 = 1;
    predVarDecisionGiven2 = 1;
end

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
function logLikelihoodCompare = PoissonLogLikelihoodCompare(theResponses,meanResponses)

logLikelihoodCompare = sum(-meanResponses + log(meanResponses).*theResponses);

end