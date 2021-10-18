function pCorrect = SupportVectorMachineObserverNAlternativeFC(meanResponses, nTestTrials)
% SVM-based Monte-Carlo simulation of probability correct for N-alternative forced choice
%
% Synopsis:
%    probCorrect = SupportVectorMachineObserverNAlternativeFC(meanResponses,nSimulatedTrials)
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
%    nSimulatedTrials -    Number of Monte-Carlo (SVM training and SVM testing) 
%                          simulated trials to use.
% Outputs:
%    probCorrect -         Probability correct
%
% Optional Key-Value Pairs

% History
%    09/21/21   NPC  Wrote it.

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
    base = 100;
    deltas = [0 1 2 5 10 20 40 80];
    nSimulatedTrials = 1000;
    for i = 1:length(deltas)
        meanResponses = [ [base+deltas(i) base base base]', [base base+deltas(i) base base]', [base base base+deltas(i) base]', [base base base base+deltas(i) ]'];
        pCorrects(i) = PoissonIdealObserverNAlternativeFC(meanResponses,nSimulatedTrials);
        pCorrects2(i) = SupportVectorMachineObserverNAlternativeFC(meanResponses,nSimulatedTrials);
        fprintf('Four alternative FC, base = %d, delta = %d, pCorrect = %0.2g\n', base, deltas(i), pCorrects(i));
    end
    figure; clf; hold on;
    plot(deltas,pCorrects,'ro-','MarkerFaceColor','r','MarkerSize',12);
    plot(deltas,pCorrects2,'ks-','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',10);
    xlabel('Delta');
    ylabel('pCorrect');
    xlim([0 80]); ylim([0.2 1]);
    set(gca, 'FontSize', 14);
    grid on; box on;
    legend({'Monte Carlo','multi-class SVM'},'Location','NorthWest');
    title(sprintf('4AFC Poisson Ideal Observer, base = %d',base));
%}
%{
    % Compare to analytic ideal observer for TAFC, where the normal
    % approximation to the decision variable is used.  Clearly doing
    % about the same thing, but not exactly the same.
    clear; close all;
    dimension = 10;
    base = 10;
    deltas = 0.04*[0 1 2 5 10 25 50 100];
    nSimulatedTrials = 2000;
    for i = 1:length(deltas)
        %meanResponses = [ [base+deltas(i) base base base]', [base base+deltas(i) base base]'];
        % Set up mean resopnses for each stimulus
        meanResponse1 = base*ones(dimension, 1);
        meanResponse2 = meanResponse1 + deltas(i)*ones(dimension,1);

        % Concatenate for TAFC where both stimuli are presented on each trial.
        meanResponses = [[meanResponse1 ; meanResponse2] , [meanResponse2 ; meanResponse1]];

        pCorrects(i) = PoissonIdealObserverNAlternativeFC(meanResponses,nSimulatedTrials, false);
        pCorrects1(i) = analyticPoissonIdealObserver(meanResponse1,meanResponse2);
        pCorrects2(i) = SupportVectorMachineObserverNAlternativeFC(meanResponses,nSimulatedTrials);
    end
    figure; clf; hold on;
    plot(deltas,pCorrects,'ro-','MarkerFaceColor','r','MarkerSize',12);
    plot(deltas,pCorrects1,'bo-','MarkerFaceColor','b','MarkerSize',10);
    plot(deltas,pCorrects2,'ks-','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',10);
    xlabel('Delta');
    ylabel('pCorrect');
    xlim([0 max(deltas)]); ylim([0.4 1]);
    set(gca, 'FontSize', 14);
    grid on; box on;
    legend({'Monte Carlo','Analytic Approx.', 'SVM'},'Location','NorthWest');
    title(sprintf('TAFC Poisson Ideal Observer & SVM Computational Observer, base = %d',base));
%}

    % Setup
    [nSensors, nAlternatives] = size(meanResponses);
    nTrainingTrials = nTestTrials;
    nSimulatedTrials = nTrainingTrials + nTestTrials;

    % Allocate memory
    whichAlternative = zeros(nSimulatedTrials,1);
    theResponses = zeros(nSimulatedTrials, nSensors);
    
    % Compute Poisson responses
    for ii = 1:nSimulatedTrials
        % Chose alternative
        whichAlternative(ii,1) = randi(nAlternatives);

        % Compute responses
        theResponses(ii,:) = poissrnd(meanResponses(:,whichAlternative(ii,1)));
    end

    % Normalize features
    [theResponses, mu, sigma] = featureNormalize(theResponses);

    % Split into train and test data sets
    xTrain = theResponses(1:nTrainingTrials,:);
    yTrain = whichAlternative(1:nTrainingTrials);
    xTest = theResponses(nTrainingTrials+1:end,:);
    yTest = whichAlternative(nTrainingTrials+1:end);

    % Choose type of SVM 
    %template = templateSVM('KernelFunction', 'gaussian');
    %template = templateSVM('KernelFunction', 'polynomial', 'PolynomialOrder', 3);
    template = templateSVM('KernelFunction', 'linear');
    codingDesign = 'binarycomplete'; % Choose from  {'allpairs', 'binarycomplete', 'denserandom', 'onevsall', 'ordinal', 'sparserandom', 'ternarycomplete'
    
    % Train the SVM
    trainedMultiClassSVM = fitcecoc(xTrain, yTrain, 'Learners', template, 'Coding', codingDesign);
    
    % Assess performance on the test data set
    nCorrect = 0;
    for k = 1:numel(yTest)
        predictedY = predict(trainedMultiClassSVM,xTest(k,:));
        if (predictedY == yTest(k))
            nCorrect = nCorrect + 1;
        end
    end
    
    % Compute pCorrect
    pCorrect = nCorrect / numel(yTest);
    pCorrect2 = (1-loss(trainedMultiClassSVM, xTest, yTest));
end

function [X, mu, sigma] = featureNormalize(X)
    mu = mean(X,1);
    sigma = std(X,0,1);
    X = bsxfun(@minus, X, mu);
    X = bsxfun(@rdivide, X, sigma);
end
