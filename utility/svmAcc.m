function [acc, w] = svmAcc(data, labels, nFolds, svmType, opts)
%% function svmAcc(dataMatrix, labels, [nFolds], [svmType], [opts])
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
%                 'svm'       - libsvm
%                 'linear'    - liblinear (default)
%    opts       - string, options to be supplied to svm packages
%
%  Outputs:
%    acc        - accuracy estimated by cross validation
%    w          - svm coefficients (linear kernel only)
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
%  See also:
%    svmtrain, svmpredict, train, predict
%
%  HJ, ISETBIO Team, 2015

%% Check Inputs
if nargin < 1, error('Data matrix required'); end
if nargin < 2, error('Labels for data required'); end
if nargin < 3, nFolds  = 5; end
if nargin < 4, svmType = 'linear'; end

%% Divide data into nFolds
%  Normalize Data
%  Normalization is important when data values are small, i.e. volts image
%  If we values are large, i.e. photon absorptions, it's fine to skip this
%  step
data = (data-repmat(min(data),[size(data, 1) 1])) ...
    ./ repmat(max(data)-min(data),[size(data,1 ) 1]);

% The data should not contain any NaN inside
% If there's any NaN, set it to zero
data(isnan(data)) = 0;
data = sparse(data);

switch svmType
    case 'linear'
        % Liblinear routine
        if notDefined('opts')
            opts = sprintf('-q -B 1 -s 2 -C -v %d', nFolds);
        end
        res = train(labels, data, opts);
        if isscalar(res)  % no -C specified
            acc = res/100;
            % the following lines are not right, should be fixed. HJ
            if nargout > 1
                s = train(labels, data, '-q -B 1 -s 2');
                w = s.w(1:end-1);
            end
        elseif length(res) == 2
            % -C is specified
            c = res(1);
            acc = res(2);
            % the following lines are not correct, should be fixed. HJ
            if nargout > 1
                opts = sprintf('-q -B 1 -s 2 -c %f', c);
                s = train(labels, data, opts);
                w = s.w(1:end-1);
            end
        end
    case 'svm'
        error('NYI');
    otherwise
        error('Unknown svmType');
end