function [acc, w] = svmAcc(data, labels, nFolds, svmType, opts)
% Compute the svm prediction accuracy
%
% Syntax:
%   [acc, w] = svmAcc(data, labels, [nFolds], [svmType], [opts])
%
% Description:
%    The function is intended to compute the svm prediciton accuracy.
%
% Inputs:
%    data    - The data matrix, an M x N matrix containing data for M
%              instances and N features
%    labels  - An M x 1 vector containing a label for each instance
%    nFolds  - (Optional) The number of folds to be used for
%              cross-validation. Default is 5.
%    svmType - (Optional) A string to indicate which svm package will be
%              used. Default is 'linear'. Options are:
%       'svm'    - libsvm
%       'linear' - liblinear (Default)
%    opts    - (Optional) Options to be supplied to the svm packages.
%              Default is '-q -B 1 -s 2 -C -v %d', where %d is nFolds
%
% Outputs:
%    acc     - The accuracy estimated by cross validation
%    w       - The svm coefficients (only provided when classifying the
%              'linear' option).
%
% Notes:
%    * [Note: XXX - For libsvm and liblinear, there are Matlab interfaces.
%      However, for most other packages, there are only C++ implementation.
%      We might write Matlab interface some time in the future. At this
%      point, we use system() command to invoke the compiled files. To run
%      in this mode successfully, you might need to compile the source code
%      on your own.]
%    * [Note: JNM - input check is primative, we should update this.]
%    * [Note: JNM - the default options are the same for both cases, so can
%      we move the notDefined step outside of the case statement so that it
%      only has to be written once?]
%    * [Note: JNM - Differences between this and svmClassifyAcc? Aside from
%      opts? And that this one was created a year later? (This one doesn't
%      have a supported 'svm' case)]
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
%    xx/xx/15  HJ   ISETBIO Team, 2015
%    12/14/17  jnm  Formatting

%  Examples:
%{
    dataMatrix = [1 2; 3 4]
	labels = ['low', 'high']
    acc = svmAcc(dataMatrix, labels, 10)
%}

%% Check Inputs
if nargin < 1, error('Data matrix required'); end
if nargin < 2, error('Labels for data required'); end
if nargin < 3, nFolds  = 5; end
if nargin < 4, svmType = 'linear'; end

%% Divide data into nFolds
% Normalize Data - Normalization is important when data values are small,
% i.e. volts image If we values are large, i.e. photon absorptions, it's
% fine to skip this step
data = (data-repmat(min(data),[size(data, 1) 1])) ...
    ./ repmat(max(data) - min(data),[size(data,1 ) 1]);

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
            acc = res / 100;
            % the following lines are not right, should be fixed. HJ
            if nargout > 1
                s = train(labels, data, '-q -B 1 -s 2');
                w = s.w(1:end - 1);
            end
        elseif length(res) == 2
            % -C is specified
            c = res(1);
            acc = res(2);
            % the following lines are not correct, should be fixed. HJ
            if nargout > 1
                opts = sprintf('-q -B 1 -s 2 -c %f', c);
                s = train(labels, data, opts);
                w = s.w(1:end - 1);
            end
        end
    case 'svm'
        error('NYI');
    otherwise
        error('Unknown svmType');
end