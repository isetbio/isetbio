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
    % Compare to analytic ideal observer for TAFC, where the normal
    % approximation to the decision variable is used.  Clearly doing
    % about the same thing, but not exactly the same.
    clear; close all;
    base = 1000;
    deltas = [0 1 2 5 10 20 40 80];
    nSimulatedTrials = 10000;
    for i = 1:length(deltas)
        meanResponses = [ [base+deltas(i) base base base]', [base base+deltas(i) base base]'];
        pCorrects(i) = PoissonIdealObserverNAlternativeFC(meanResponses,nSimulatedTrials);
        pCorrects1(i) = analyticPoissonIdealObserver(meanResponses(:,1),meanResponses(:,2));
        pCorrects2(i) = SupportVectorMachineObserverNAlternativeFC(meanResponses,nSimulatedTrials);
    end
    figure; clf; hold on;
    plot(deltas,pCorrects,'ro-','MarkerFaceColor','r','MarkerSize',12);
    plot(deltas,pCorrects1,'bo-','MarkerFaceColor','b','MarkerSize',10);
    plot(deltas,pCorrects2,'ks-','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',10);
    xlabel('Delta');
    ylabel('pCorrect');
    xlim([0 80]); ylim([0.4 1]);
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
