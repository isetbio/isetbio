function acc = svmClassifyAcc(dataMatrix, labels, nFolds, svmType, opts)
%% function getSVMAccuracy(dataMatrix, labels, [nFolds], [svmType], [opts])
%    compute the svm prediction accuracy
%
%  Supported packages:
%    libsvm: free copy can be downloaded from
%            http://www.csie.ntu.edu.tw/~cjlin/libsvm/
%    liblinear: free copy can be downloaded from
%            http://www.csie.ntu.edu.tw/~cjlin/liblinear/
%
%  Inputs:
%    dataMatrix - M-by-N matrix, containing data for M instances and N
%                 features
%    labels     - M-by-1 vector, label for each instance
%    nFolds     - number of folds for cross validation
%    svmType    - string, indicate which svm package to be used:
%                 'svm'       - libsvm (default)
%                 'linear'    - liblinear (linear kernal)
%                 'ranksvm'   - liblinear - ranksvm (fast / scalable mode)
%                 'warmstart' - liblinear - warmstart (incremental mode)
%    opts       - string, options to be supplied to svm packages
%
%  Outputs      - 2-by-1 vector, containing average accuray and standard
%                 deviation
%
%  Example:
%    acc = svmClassifyAcc(dataMatrix, labels, 10, 'linear')
%
%  Notes:
%    For libsvm and liblinear, there are Matlab interfaces. However, for
%    most other packages, there are only C++ implementation. We might write
%    Matlab interface some time in the future. At this point, we use
%    system() command to invoke the compiled files. To run in this mode
%    successfully, you might need to compile the source code on your own
%
%  TODO:
%    There are other optimized packages for non-linear kernal svm and we
%    might add them some time in the future. Some packages are listed as
%    below:
%      1. Core Vector Machine
%      2. FaLKM-lib (fast, local kernal svm)
%
%  See also:
%    svmtrain, svmpredict, train, predict
%
%  (HJ) ISETBIO Team, 2014

%% Check Inputs
if nargin < 1, error('Data matrix required'); end
if nargin < 2, error('Labels for data required'); end
if nargin < 3, nFolds  = 5; end
if nargin < 4, svmType = 'svm'; end

%% Divide data into nFolds
%  Normalize Data
%  Normalization is important when data values are small, i.e. volts image
%  If we values are large, i.e. photon absorptions, it's fine to skip this
%  step
dataMatrix = (dataMatrix-repmat(min(dataMatrix),[size(dataMatrix, 1) 1])) ...
    ./ repmat(max(dataMatrix)-min(dataMatrix),[size(dataMatrix,1 ) 1]);

% The data should not contain any NaN inside
% If there's any NaN, set it to zero
dataMatrix(isnan(dataMatrix)) = 0;

%  Random permute the data
[M, ~] = size(dataMatrix);
ind   = randperm(M);

% Train and Test
accHistory  = zeros(nFolds,1);
instPerFold = round(M/nFolds);
for i = 1 : nFolds
    if i < nFolds
        trainIndx = [ind(1:(i-1)*instPerFold) ...
            ind(i*instPerFold+1:end)];
        testIndx  = ind((i-1)*instPerFold+1:i*instPerFold);
    else
        trainIndx = ind(1:(i-1)*instPerFold);
        testIndx  = ind((i-1)*instPerFold+1:end);
    end
    trainData = sparse(dataMatrix(trainIndx,:));
    testData  = sparse(dataMatrix(testIndx,:));
    % Train
    switch svmType
        case 'linear'
            % Liblinear routine
            if notDefined('opts'), opts = '-s 2 -q'; end
            svmStruct = train(labels(trainIndx),trainData,opts);
            [~,curAcc,~] = predict(labels(testIndx),testData,svmStruct,'-q');
        case 'svm'
            % LibSVM routine
            % Parameters explaination:
            %   http://www.csie.ntu.edu.tw/~cjlin/libsvm/
            if notDefined('opts'), opts = '-s 2 -q'; end
            svmStruct = svmtrain(labels(trainIndx),trainData,opts);
            [~,curAcc,~] = svmpredict(labels(testIndx), testData, ...
                svmStruct,'-q');
        case 'ranksvm'
            if notDefined('opts'), opts = '-s 2 -q'; end
            libsvmwrite('ranksvm.train', labels(trainIndx), trainData);
            libsvmwrite('ranksvm.test', labels(testIndx), testData);
            cmd = fullfile(ranksvmRootPath, 'train');
            cmd = sprintf('%s %s %s model', cmd, opts, 'ranksvm.train');
            system(cmd);
            
            % now model is saved to ranksvm.train.model
            % do prediction
            cmd = fullfile(ranksvmRootPath, 'predict');
            cmd = sprintf('%s -q %s %s %s', cmd, 'ranksvm.test', ...
                            'model', 'pred_labels');
            system(cmd);
            
            % load back
            pred_labels = importdata('pred_labels');
            curAcc = sum(pred_labels == labels(testIndx)) / ...
                            length(pred_labels) * 100;
                        
            % clean up
            delete('ranksvm.train', 'ranksvm.test');
            delete('model', 'pred_labels');
        case 'increment'
        otherwise
            error('Unknown svm type');
    end
    accHistory(i) = curAcc(1) / 100; % Convert to between 0~1
end

% Report average and std
acc(1) = mean(accHistory);
acc(2) = std(accHistory);

end