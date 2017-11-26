function wvfParams = wvfComputeOptimizedConePSF(wvfParams)
% Optimize the PSF by the cones with sensitivities, weighting, & criterion
%
% Syntax:
%   wvfParams = wvfComputeOptimizedConePSF(wvfParams)
%
% Description:
%    Optimize the PSF seen by the cones, given the cone sensitivities, a
%    weighting spectral power distribution, and a criterion. Optimization
%    is performed on the defocus parameter. 
%
% Inputs:
%    wvfParams - 
%
% Outputs:
%    wvfParams - 
%
% Notes:
%    * [NOTE: DHB - This function is under development. The idea is that
%       the amound of defocus that produces the best PSF, as seen by a
%       particular cone class, is not determined trivially and is best
%       found via numerical optimization. We may never need this, and
%       certainly don't need it right now. I am moving to an
%       "underDevelopment_wavefront" directory and putting in an error
%       message at the top so that people don't think it might work.]
%    * [NOTE: JNM - The example is broken, but at least 'semi-present' now]
%

% History:
%    08/26/11  dhb  Wrote it.
%    08/29/11  dhb  Don't need to center or circularly average here.
%              dhb  Print warning if optimal value is at search bound.
%    09/07/11  dhb  Rename. Use wvfParams for i/o.
%	 11/14/17  jnm  Comments & formatting

% Examples:
%{
    %% Load cone sensitivities, set weighting spectrum.
    S = [400 5 61];
    wls = SToWls(S);
    load T_cones_ss2;
    T_cones = SplineCmf(S_cones_ss2,T_cones_ss2,S);
    load spd_D65
    weightingSpectrum = SplineSpd(S_D65,spd_D65,S);

    % Specify datafile for Zernike coefficients
    zernikeFile = 'sampleZernikeCoeffs.txt';
    measpupilMM = 6;
    theZernikeCoeffs = importdata(zernikeFile);

    wvfParams0 = wvfCreate;
    wvfParams0 = wvfSet(wvfParams0,'measured pupil',measpupilMM);
    wvfParams0 = wvfSet(wvfParams0,'calculated pupil',3);
    wvfParams0 = wvfSet(wvfParams0,'wave',wls);
    wvfParams0 = wvfSet(wvfParams0,'defocusDiopters',0);
    wvfParams0 = wvfSet(wvfParams0,'fieldSampleSize',16.212/201);
    wvfParams0 = wvfSet(wvfParams0,'fieldsizemm',16.212);
    wvfParams0.T_cones = T_cones;
    whichRow = floor(wvfGet(wvfParams0,'npixels')/2) + 1;
    wvfParams0 = wvfSet(wvfParams0,'sceparams',sceCreate(wls,'none'));
    wvfParams0 = wvfComputePSF(wvfParams0);
    wvfParams0.coneWeights = [1 1 0];
    wvfParams0.criterionFraction = 0.9;

    wvfParams0 = wvfComputeOptimizedConePSF(wvfParams0)
%}

%% This is not yet working.
error('This function is under development and not yet working');

options = optimset('fmincon');
options = optimset(options, 'Diagnostics', 'off', 'Display', 'off', ...
    'LargeScale', 'off', 'Algorithm', 'active-set');
%options = optimset(options, 'TypicalX', 0.1, 'DiffMinChange', 1e-3);

% If the parallel toolbox is present, check and use if available
if exist('matlabpool', 'builtin')
    if (exist('IsCluster', 'file') && IsCluster && matlabpool('size') > 1)
        options = optimset(options, 'UseParallel', 'always');
    end
end

% Initial defocus and bounds (diopters)
diopterBound = 4;
x0 = 0;
vlb = -diopterBound;
vub = -vlb;

% Optimize focus
x = fmincon(@InlineMinFunction, x0, [], [], [], [], vlb, vub, [], options);

% Set up return values
defocusDiopters = x;
if (abs(defocusDiopters) >= diopterBound)
    fprintf(['WARNING: defocus found in wvfComputeOptimizedConePSF is '...
        'at search limit of %0.1f diopters\n'], diopterBound)
end
[f, tmpWvfParams] = InlineMinFunction(defocusDiopters);
wvfParams = tmpWvfParams;

    function [f, tmpWvfParams] = InlineMinFunction(x)
        % Calculate the inline minumum
        %
        % Syntax:
        %   [f, tmpWvfParams] = InlineMinFunction(x)
        %
        % Description:
        %    This private function calculates the inline minimum and
        %    returns the wavefront struct and the focus.
        %
        % Inputs:
        %    x            - Wavefront struct defocus diopters
        %
        % Outputs:
        %    f            - Focus?
        %    tmpWvfParams - wavefront struct temporary variable.
        %

        tmpWvfParams = wvfParams;
        tmpWvfParams.defocusDiopters = x;
        tmpWvfParams = wvfComputeConePSF(tmpWvfParams);
        nCones = size(tmpWvfParams.T_cones, 1);
        f = 0;
        for j = 1:nCones
            %temppsf = psfCircularlyAverage(psfCenter(conepsf(:, :, j)));
            critRadius(j) = psfFindCriterionRadius(...
                tmpWvfParams.conepsf(:, :, j), ...
                tmpWvfParams.criterionFraction);
            f = f + tmpWvfParams.coneWeights(j) * critRadius(j);
        end
        
    end
end
