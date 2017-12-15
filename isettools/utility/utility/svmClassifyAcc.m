function [acc, w] = svmClassifyAcc(data, labels, nFolds, svmType, opts)
% Classify the svm prediction accuracy
%
% Syntax:
%   [acc, w] = svmClassifyAcc(data, labels, [nFolds], [svmType], [opts])
%
% Description:
%    The function is intended to classify the svm prediciton accuracy.
%
% Inputs:
%    data    - The data matrix, an M x N matrix containing data for M
%              instances and N features
%    labels  - An M x 1 vector containing a label for each instance
%    nFolds  - (Optional) The number of folds to be used for
%              cross-validation. Default is 5.
%    svmType - (Optional) A string to indicate which svm package will be
%              used. Default is 'svm'. Options are:
%       'svm'    - libsvm (Default)
%       'linear' - liblinear (linear kernel)?
%    opts    - (Optional) Options to be supplied to the svm packages.
%              Default is '-s 2 -q'
%
% Outputs:
%    acc     - A 2 x 1 vector, containing the average accuracy and the
%              standard deviation.
%    w       - The svm coefficients (only provided when classifying the
%              'linear' option).
%
% Notes:
%    * [Note: JNM - input check is primative, we should update this.]
%    * [Note: JNM - svmtrain is about to be deprecated, we need to update
%      the use of this function or this will break shortly.]
%    * [Note: JNM - the default options are the same for both cases, so can
%      we move the notDefined step outside of the case statement so that it
%      only has to be written once?]
%    * [Note: JNM - Differences between this and svmAcc? Aside from opts?]
%
% References: 
%    Free copies of the following supported packages can be downloaded at:
%       libsvm: http://www.csie.ntu.edu.tw/~cjlin/libsvm/
%       liblinear: http://www.csie.ntu.edu.tw/~cjlin/liblinear/
%
%  See Also:
%    svmtrain, svmpredict, train, predict
%

% History:
%    xx/xx/14  HJ   ISETBIO Team, 2014
%    12/14/17  jnm  Formatting

%  Examples:
%{
    dataMatrix = [1 2; 3 4]
	labels = ['low', 'high']
    acc = svmClassifyAcc(dataMatrix, labels, 10, 'linear')
%}


%% Check Inputs
warning('The cross-validation part is out of date. Should use -v opts');
if nargin < 1, error('Data matrix required'); end
if nargin < 2, error('Labels for data required'); end
if nargin < 3, nFolds = 5; end
if nargin < 4, svmType = 'svm'; end

%% Divide data into nFolds
% Normalize Data - Normalization is important when data values are small, 
% i.e. volts image If we values are large, i.e. photon absorptions, it's
% fine to skip this step
data = (data - repmat(min(data), [size(data, 1) 1])) ...
    ./ repmat(max(data) - min(data), [size(data, 1 ) 1]);

% zero mean
data = 2 * bsxfun(@minus, data, mean(data));

% The data should not contain any NaN inside
% If there's any NaN, set it to zero
data(isnan(data)) = 0;

%  Random permute the data
[M, ~] = size(data);
ind = randperm(M);

% Train and Test
accHistory = zeros(nFolds, 1);
instPerFold = round(M / nFolds);
w = zeros(size(data, 2), nFolds);
for i = 1 : nFolds
    if i < nFolds
        trainIndx = [ind(1:(i - 1) * instPerFold) ...
            ind(i * instPerFold + 1:end)];
        testIndx = ind((i - 1) * instPerFold + 1:i * instPerFold);
    else
        trainIndx = ind(1:(i - 1) * instPerFold);
        testIndx = ind((i - 1) * instPerFold + 1:end);
    end
    trainData = sparse(data(trainIndx, :));
    testData = sparse(data(testIndx, :));
    % Train
    switch svmType
        case 'linear'
            % Liblinear routine
            if notDefined('opts'), opts = '-s 2 -q'; end
            svmStruct = train(labels(trainIndx), trainData, opts);
            [~, curAcc, ~] = predict(labels(testIndx), testData, ...
                svmStruct, '-q');
            w(:, i) = svmStruct.w;
        case 'svm'
            % LibSVM routine
            % Parameters explaination:
            %   http://www.csie.ntu.edu.tw/~cjlin/libsvm/
            if notDefined('opts'), opts = '-s 2 -q'; end
            svmStruct = svmtrain(labels(trainIndx), trainData, opts);
            [~, curAcc, ~] = svmpredict(labels(testIndx), testData, ...
                svmStruct, '-q');
        otherwise
            error('Unknown svm type');
    end
    accHistory(i) = curAcc(1) / 100; % Convert to between 0~1
end

% Report average and std
acc(1) = mean(accHistory);
acc(2) = std(accHistory);

end