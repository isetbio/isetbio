function [primary,predictedSpd,errorFraction,gamutMargin] = SpdToPrimary(calOrCalStruct, targetSpd, varargin)
% Converts a spectrum into normalized primary OneLight mirror values.
%
% Syntax:
%     primary = SpdToPrimary(calOrCalStruct, targetSpd)
%     primary = SpdToPrimary(calOrCalStruct, targetSpd, 'lambda', 0.01)
%     primary = SpdToPrimary(calOrCalStruct, targetSpd, 'verbose, true)
%
% Description:
%    Convert a spectral power distribution to the linear 0-1 fraction of
%    light that we need from each primary.
%
%    No gamma correction is applied to the primary values.
%
%    This routine also allows for a 'differentialMode' which is false unless
%    the 'differentialMode' key value pair is passed.
%
%    This routine keeps values in the range [0,1] in normal mode, and in
%    range [-1,1] in differential mode.
%
%    The routine works by using lsqlin to minimize the SSE between target and
%    desired spectra.  The value of the 'lambda' key is smoothing parameter.
%    This weights an additional error term that tries to minimize the SSE
%    of the difference between neighboring primary values. This can reduce
%    ringing in the obtained primaries, at the cost of increasing the SSE
%    to which the target spd is reproduced.
%
%    Set value of 'checkSpd' to true to force a check on how well the
%    target is acheived.
%
%    This takes vector, not matrix, input. Put a loop around the call to
%    this function, or edit so that it loops over the columns of a passed
%    matrix, if you want to apply to more than on spd.
%
%    Also note that if you know that the target is a linear combination of
%    the device primaries, you can do this much faster with:
%       calOrCalStruct = ObjectToHandleCalOrCalStruct(calOrCalStruct);
%       P_device = calOrCalStruct.get('P_device');
%       P_ambient = calOrCalStruct.get('P_ambient');
%       primary = P_device\(targetSpd-P_ambient).
%
% Inputs:
%    calOrCalStruct    - Calibration struct or object.
%    targetSpd         - Column vector providing the target spectrum, sampled at the wavelengths
%                        used in the calibration file (typically 380 to 780 nm in 2 nm steps)..
%
% Outputs:
%    primary           - Column vector containing the primary values for each effective primary
%                        of the OneLight. nPrimaries is the number of
%                        effective primaries. Not gamma corrected.
%    predictedSpd      - The spd predicted for the returned primaries.
%    errorFraction     - How close the predictedSpd came to the target, in
%                        fractional terms.
%    gamutMargin       - How far out of gamut are primaries. Positive means
%                        out of gamut, negative is how far in gamut.
%
% 
% Optional Key-Value Pairs:
%  'verbose'           - Boolean (default false). Provide more diagnostic output.
%  'lambda'            - Scalar (default 0.005). Value of primary smoothing
%                        parameter.  Smaller values lead to less smoothing,
%                        with 0 doing no smoothing at all.
%   'primaryHeadroom'  - Scalar.  Headroom to leave on primaries.  Default
%                        0. How much headroom to protect in definition of
%                        in gamut.  Range used for check and truncation is
%                        [primaryHeadroom 1-primaryHeadroom]. Do not change
%                        this default.  Sometimes assumed to be true by a
%                        caller.
%   'primaryTolerance  - Scalar. Truncate to range [0,1] if primaries are
%                        within this tolerance of [0,1]. Default 1e-6, and
%                        'checkPrimaryOutOfRange' value is true.
%   'checkPrimaryOutOfRange' - Boolean. Perform primary tolerance check.
%                        Default true.
%   'differentialMode' - Boolean. Run in differential
%                       mode.  This means, don't subtract dark light.
%                       Useful when we want to find delta primaries that
%                       produce a predicted delta spd. Default false.
%   'checkSpd'        - Boolean (default false). Because of smoothing and
%                       gamut limitations, this is not guaranteed to
%                       produce primaries that lead to the predictedSpd
%                       matching the targetSpd.  Set this to true to check
%                       force an error if difference exceeds tolerance.
%                       Otherwise, the toleranceFraction actually obtained
%                       is retruned. Tolerance is given by spdTolerance.
%   'spdToleranceFraction' - Scalar (default 0.01). If checkSpd is true, the
%                       tolerance to avoid an error message is this
%                       fraction times the maximum of targetSpd.
%   'whichSpdToPrimaryMin' - String, what to minimize (default 'leastSquares')
%                           * 'leastSquares' Mimimize sum of squared error,
%                             respecting lambda parameter as well.  Fast.
%                           * 'fractionalError' Minimize fractional squared
%                              error. Way slower than 'leastSquares'.  This
%                              method also respects lambda constraint, but
%                              the relative scale of the error is different
%                              so you need to adjust lambda by hand. This
%                              exists mainly for debugging. Not recommended
%                              for everyday use.
%   'maxSearchIter'   - Control how long the search goes for, when using
%                       'fractionalError' method. Default, 300.  Reduce if you
%                       don't need to go that long and things will get
%                       faster.
%
% See also:
%   PrimaryToSpd, PrimaryToSettings
%

% History:
%    09/11/21  dhb  Wrote from OLSpdToPrimary.
 
% Examples:
%{
    % ETTBSkip
    %
    % This example needs updating.  Still written for ancestral routine
    % OLSpdToPrimary.
    cal = OLGetCalibrationStructure('CalibrationType','DemoCal','CalibrationFolder',fullfile(tbLocateToolbox('OneLightToolbox'),'OLDemoCal'),'CalibrationDate','latest');
    primaryIn = rand(size(cal.computed.pr650M,2),1);
    spd1 = OLPrimaryToSpd(cal,primaryIn);
    primaryOut = OLSpdToPrimary(cal,spd1,'primaryHeadroom',0,'lambda',0);
    spd2 = OLPrimaryToSpd(cal,primaryOut);
    figure; clf;
    plot(primaryIn,primaryOut,'ro','MarkerSize',4,'MarkerFaceColor','r');
    axis([0 1 0 1]); axis('square');
    figure; clf; hold on;
    plot(spd1,'r','LineWidth',4);
    plot(spd2,'b','LineWidth',2);
%}

%% Parse the input
p = inputParser;
p.addParameter('verbose', false, @islogical);
p.addParameter('lambda', 0.005, @isscalar);
p.addParameter('primaryHeadroom', 0, @isscalar);
p.addParameter('primaryTolerance', 1e-6, @isscalar);
p.addParameter('checkPrimaryOutOfRange', true, @islogical);
p.addParameter('differentialMode', false, @islogical);
p.addParameter('checkSpd', false, @islogical);
p.addParameter('whichSpdToPrimaryMin', 'leastSquares', @ischar);
p.addParameter('spdToleranceFraction', 0.01, @isscalar);
p.addParameter('maxSearchIter',300,@isscalar);
p.parse(varargin{:});
params = p.Results;

%% Make sure we have @CalStruct object that will handle all access to the calibration data.
%
% From this point onward, all access to the calibration data is accomplised via the calStructOBJ.
[calStructOBJ, inputArgIsACalStructOBJ] = ObjectToHandleCalOrCalStruct(calOrCalStruct);
if (~inputArgIsACalStructOBJ)
    error('The input (calOrCalStruct) is not a cal struct.');
end

%% Parameters
if (p.Results.verbose)
    fminconDisplaySetting = 'iter';
else
    fminconDisplaySetting = 'off';
end

%% Check wavelength sampling
S = calStructOBJ.get('S');
wls = SToWls(S);
nWls = size(targetSpd,1);
if (nWls ~= S(3))
    error('Wavelength sampling inconsistency between passed spectrum and calibration');
end

%% Dark light
% 
% In differential mode, we ignore the dark light, otherwise we snag it from
% the calibration file.
if params.differentialMode
    darkSpd = zeros(size(calStructOBJ.get('P_ambient')));
else
    darkSpd = calStructOBJ.get('P_ambient');
end

%% Device primaries
P_device = calStructOBJ.get('P_device');

%% Find primaries the linear way, without any constraints
%
% This would be the most straightforward way to find the primaries, but for many
% spectra the primary values really ring, which is why we use the search based
% method below and enforce a smoothing regularization constraint.
%
% We skip this step unless we are debuging.
DEBUG = 0;
if (DEBUG)
    pinvprimary = pinv(P_device) * (targetSpd - darkSpd);
    if params.verbose
        fprintf('Pinv values: min = %g, max = %g\n', min(pinvprimary(:)), max(pinvprimary(:)));
    end
end

%% Initialize some primaries for search
initialPrimary = 0.5*ones(size(P_device,2),1);

%% Linear constraint setup
%
% This first constraint (C1,d1) minimizes the error between the predicted spectrum
% and the desired spectrum.
C1 = P_device;
d1 = targetSpd - darkSpd;

% Scale so that error is in a known range. Prevents lsqlin from
% thinking it has a good solution when units of light lead to small
% numbers in the target spectrum.
C1 = C1;
d1 = d1;

% The second constraint computes the difference between between neighboring
% values and tries to make this small.  How much this is weighted 
% depends on the value of params.lambda.  The bigger params.lambda, the
% more this constraint kicks in.
nPrimaries = size(P_device,2);
C2 = zeros(nPrimaries -1, nPrimaries );
for i = 1:nPrimaries -1
    C2(i,i) = params.lambda;
    C2(i,i+1) = -params.lambda;
end
d2 = zeros(nPrimaries-1,1);

% Paste together the target and smoothness constraints
%
% Scale so that SSE for target is in reasonable range. 
% Scale lambda term as well, to preserve meaning of lambda
% vis-a-vis the time before we did this scaling.
%
% Check for very small means, which can happen in differential mode,
% and put in a reasonable numerical value. This value was obtained from
% the dark measurement level in the demo calibration file.
meanScale = abs(mean(targetSpd));
if (meanScale < 5e-5)
    meanScale = 5e-5;
end
C = [C1 ; C2]/meanScale;
d = [d1 ; d2]/meanScale;

%% Primary bounds for searches
if params.differentialMode
    vub = ones(size(initialPrimary))  - p.Results.primaryHeadroom;
    vlb = -1*ones(size(initialPrimary)) + p.Results.primaryHeadroom;
else
    vub = ones(size(initialPrimary))  - p.Results.primaryHeadroom;
    vlb = zeros(size(initialPrimary)) + p.Results.primaryHeadroom;
end

%% Search for primaries that hit target as closely as possible,
% 
% Subject to lambda weighting.  This is slower and less accurate than
% lsqlin, presumably because lsqlin is written to know about the linear,
% structure of the error.
%
% But, fmincon gives more control so keeping this here in case we ever want 
% explore further;

% options = optimset('fmincon');
% options = optimset(options,'Diagnostics','off','Display',fminconDisplaySetting,'LargeScale','off','Algorithm','active-set', 'MaxIter', p.Results.maxSearchIter, 'MaxFunEvals', 100000, 'TolFun', 1e-7, 'TolCon', 1e-6, 'TolX', 1e-4);
% fminconx = fmincon(@(x) OLFindSpdFun(x, cal, targetSpd, p.Results.lambda, p.Results.differentialMode), ... 
%     initialPrimary,[],[],[],[],vlb,vub, ...
%     [], ...
%     options);
% fminconResnorm = FindSpdFun(fminconx, cal, targetSpd, p.Results.lambda, p.Results.differentialMode);

switch (p.Results.whichSpdToPrimaryMin)
    case 'fractionalError'
        % Search minimizing fractional error rather than sum of squared error
        options = optimset('fmincon');
        options = optimset(options,'Diagnostics','off','Display',fminconDisplaySetting,'LargeScale','off','Algorithm','active-set', 'MaxIter', p.Results.maxSearchIter, 'MaxFunEvals', 100000, 'TolFun', p.Results.spdToleranceFraction/10, 'TolCon', 1e-6, 'TolX', 1e-6);
        fminconfractionalx = fmincon(@(x) FindSpdFunFractional(x, calStructOBJ, targetSpd, p.Results.lambda, p.Results.differentialMode), ...
            initialPrimary,[],[],[],[],vlb,vub, ...
            [], ...
            options);
        primary = fminconfractionalx;
        
    case 'leastSquares'
        % Use lsqlin to find primaries.
        %
        % Call into lsqlin
        lsqoptions = optimset('lsqlin');
        lsqoptions = optimset(lsqoptions,'Diagnostics','off','Display','off');
        [lsqx,lsqResnorm] = lsqlin(C,d,[],[],[],[],vlb,vub,[],lsqoptions);
        primary = lsqx;
        
    otherwise
        error('Unknown value for ''whichMin'' key');
end

%% Compare ways of getting spectrum
% 
% Need to run alternative ways by hand to get all the lines.
%{
figure; clf; hold on
plot(targetSpd,'r','LineWidth',4);
plot(OLPrimaryToSpd(cal,lsqx,'skipAllChecks',true),'g','LineWidth',1);
plot(OLPrimaryToSpd(cal,fminconfractionalx,'skipAllChecks',true),'b','LineWidth',1);
plot(OLPrimaryToSpd(cal,fminconx,'skipAllChecks',true),'k','LineWidth',1);
%}

%% Report
if params.verbose
    fprintf('Primaries: min = %g, max = %g\n', min(primary(:)), max(primary(:)));
end

%% Make sure we enforce bounds, in case lsqlin has a bit of numerical slop
[primary,~,gamutMargin] = CheckPrimaryGamut(primary, ...
    'primaryHeadroom',p.Results.primaryHeadroom, ...
    'primaryTolerance',p.Results.primaryTolerance, ...
    'checkPrimaryOutOfRange',p.Results.checkPrimaryOutOfRange, ...
    'differentialMode',p.Results.differentialMode);

%% Predict spd, and check if specified
predictedSpd = PrimaryToSpd(calStructOBJ,primary, ...
    'differentialMode',params.differentialMode);

[~,errorFraction] = CheckSpdTolerance(targetSpd,predictedSpd, ...
    'checkSpd',p.Results.checkSpd,'spdToleranceFraction',p.Results.spdToleranceFraction);
end

function f = OLFindSpdFun(primary, calStructOBJ, targetSpd, lambda, differentialMode)

% Get the prediction.  Constraint checking is done in the constraint
% function, skipped here
predictedSpd = OLPrimaryToSpd(calStructOBJ, primary, ...
    'differentialMode', differentialMode, ...
    'skipAllChecks',true);

% Fit to desired sum of squares
f1 = sum((predictedSpd-targetSpd).^2);

% Primary smoothness penalty
primaryDiffs = diff(primary).^2;
f2 = lambda*sum(primaryDiffs);

% Final error.
f = (f1 + f2)/(mean(targetSpd)^2);

end


function f = FindSpdFunFractional(primary, cal, targetSpd, lambda, differentialMode)

% Get the prediction.  Constraint checking is done in the constraint
% function, skipped here
predictedSpd = PrimaryToSpd(cal, primary, ...
    'differentialMode', differentialMode, ...
    'skipAllChecks',true);

[~, errorFraction] = CheckSpdTolerance(targetSpd,predictedSpd, ...
    'checkSpd', false);

f1 = errorFraction;

% Primary smoothness penalty
primaryDiffs = diff(primary).^2;
f2 = lambda*sum(primaryDiffs);

% Final error.
f = (f1 + f2);

end
