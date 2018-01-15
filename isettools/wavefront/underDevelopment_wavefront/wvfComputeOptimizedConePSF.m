function wvfOut = wvfComputeOptimizedConePSF(wvfIn)
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
%    wvfIn - 
%
% Outputs:
%    wvfOut - 
%
% Notes:
%    * [NOTE: DHB - This function is under development. The idea is that
%       the amount of defocus that produces the best PSF, as seen by a
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

% From Jenn's notes when she worked on this in November 2017.
% wvfComputeOptimizedConePSF
%     - Note: Was the input/output parameter intended to be non-private?
%         I think I found the problem. wvfParams is not passed into the function directly, but is called regardless.
%     - Note: How should the inlineMinFunction be addressed?
%     - Note: Example is not working! (Trying, but still struggling)
%     - Please check that all of my commentary inside the inline function is accurate. I am only like 50% certain of most of it. I have also neglected to include an example

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

options = optimset('fmincon');
options = optimset(options, 'Diagnostics', 'off', 'Display', 'off', ...
    'LargeScale', 'off', 'Algorithm', 'active-set');
%options = optimset(options, 'TypicalX', 0.1, 'DiffMinChange', 1e-3);

% Initial defocus and bounds (diopters)
diopterBound = 4;
defocusStart = 0;
vlb = -diopterBound;
vub = -vlb;

% Optimize focus
defocusFound = fmincon(@(defocus) InlineMinFunction(defocus,wvfIn), defocusStart, [], [], [], [], vlb, vub, [], options);
if (abs(x) >= diopterBound)
    fprintf(['WARNING: defocus found in wvfComputeOptimizedConePSF is '...
        'at search limit of %0.1f diopters\n'], diopterBound)
end
[~, wvfOut] = InlineMinFunction(x,wvfIn);


    end
