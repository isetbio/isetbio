function [oiEnsemble, psfEnsemble, zCoeffs] = oiEnsembleGenerate(obj, oiSamplingGridDegs, varargin)
% Create an ensemble of optical images and psfs 
%
% Brief description
%  The varargin parameters specify the dataset and other processing
%  parameters.
%
% Input:
%   obj:  A cMosaic
%
% Optional key/val pairs
%  zernikeDataBase           - Either one of {'Polans2015', 'Artal2012'}
%  subjectID                 - Subject number in the data base. Default 6.
%  pupilDiameterMM           - Pupil diameter in millimeters (default 3mm)
%  wavefrontSpatialSamples   - Default is 301 samples
%  subtractCentralRefraction - Fixes up the data a bit. Default false.
%  zeroCenterPSF             - Default true
%  deNoisedZernikeCoefficients - Seems deprecated to me (BW).
%  flipPSFUpsideDown         - Default true (aligns with image)
%  defocusMicrons            - Add this value to defocus zcoeff
%
% Outputs:
%   oiEnsemble   -  Cell array of OIs
%   psfEnsemble  -  Cell array of PSFs
%   zCoeffs      -  Zernike polynomial coeffcieints
%
% Description:
%   In which we explain more about the processing parameters.  NC and BW to
%   do together.
%
%
% See also
%   cMosaic (main class)
%

% Help
if (ischar(oiSamplingGridDegs)) && (strcmp(oiSamplingGridDegs,'help'))
    doc('cMosaic.oiEnsembleGenerate');
    return;
end

% Parse input
p = inputParser;
p.addRequired('obj', @(x)(isa(x,'cMosaic')));
p.addRequired('oiSamplingGridDegs', @(x)(isnumeric(x) && (size(x,2) == 2) && (all(isreal(x)))));
p.addParameter('zernikeDataBase', 'Polans2015', @(x)(ismember(x, {'Polans2015', 'Artal2012'})));
p.addParameter('warningInsteadOfErrorForBadZernikeCoeffs', false, @islogical);
p.addParameter('subjectID', 6, @isscalar);
p.addParameter('pupilDiameterMM', 3.0, @isscalar);
p.addParameter('wavefrontSpatialSamples', 301, @isscalar);
p.addParameter('subtractCentralRefraction', false, @islogical);
p.addParameter('zeroCenterPSF', true, @islogical);
p.addParameter('flipPSFUpsideDown', true, @islogical);
p.addParameter('upsampleFactor', [], @(x)(isempty(x) || ((isnumeric(x))&&(numel(x)==1)&&(x>0))));
p.addParameter('refractiveErrorDiopters', 0, @isnumeric);
p.parse(obj, oiSamplingGridDegs, varargin{:});

oiSamplingGridDegs = p.Results.oiSamplingGridDegs;
zernikeDataBase = p.Results.zernikeDataBase;
pupilDiamMM = p.Results.pupilDiameterMM;
subjectID = p.Results.subjectID;
subtractCentralRefraction = p.Results.subtractCentralRefraction;
wavefrontSpatialSamples = p.Results.wavefrontSpatialSamples;
zeroCenterPSF = p.Results.zeroCenterPSF;
flipPSFUpsideDown = p.Results.flipPSFUpsideDown;
upSampleFactor = p.Results.upsampleFactor;
warningInsteadOfErrorForBadZernikeCoeffs = p.Results.warningInsteadOfErrorForBadZernikeCoeffs;

% Generate the oiEnsemble
oiNum = size(oiSamplingGridDegs,1);
oiEnsemble = cell(1, oiNum);
psfEnsemble = cell(1, oiNum);

switch (zernikeDataBase)
    
    case 'Artal2012'
        % Make sure refractive error is zero, because Artal version of
        % oiForSubjectAtEccentricity doesn't understand the
        % 'refractiveErrorMicrons' key/value pair.
        if (p.Results.refractiveErrorDiopters ~= 0)
            error('Artal optics does not currently accept refractiveErrorDiopters key/value pair');
        end

        % Artal optics
        for oiIndex = 1:oiNum
            %fprintf('Generating %s optics for eccentricity: %2.1f,%2.1f degs (um/deg):%2.1f\n', ...
            %    zernikeDataBase, oiSamplingGridDegs(oiIndex,1), oiSamplingGridDegs(oiIndex,2), obj.micronsPerDegree);
            targetEcc = oiSamplingGridDegs(oiIndex,:);
            
            
            if (targetEcc(2) ~= 0)
                fprintf(2,'Artal optics not available off the horizontal meridian. Computing for vEcc = 0\n');
                targetEcc(2) = 0;
            end
            
            [theOI, thePSF, psfSupportMinutesX, psfSupportMinutesY, psfSupportWavelength, zCoeffs] = ArtalOptics.oiForSubjectAtEccentricity(subjectID, ...
                obj.whichEye, targetEcc, pupilDiamMM, obj.wave, obj.micronsPerDegree, ...
                'wavefrontSpatialSamples', wavefrontSpatialSamples, ...
                'subtractCentralRefraction', subtractCentralRefraction, ...
                'zeroCenterPSF', zeroCenterPSF, ...
                'flipPSFUpsideDown', flipPSFUpsideDown, ...
                'upsampleFactor', upSampleFactor);
            
            if (isempty(theOI))
                if (warningInsteadOfErrorForBadZernikeCoeffs)
                    fprintf(2,'Bad Zernike coefficents for the %s of Artal subject %d. Choose another subject/eye\n', obj.whichEye, subjectID);
              
                    oiEnsemble = [];
                    psfEnsemble = []; 
                    zCoeffs = [];
                    return;
                else
                    error('Bad Zernike coefficents for the %s of Artal subject %d. Choose another subject/eye', obj.whichEye, subjectID);
                end
                
            end
            
            oiEnsemble{oiIndex} = theOI;
            psfEnsemble{oiIndex} = struct(...
                'data', thePSF, ...
                'supportX', psfSupportMinutesX, ...
                'supportY', psfSupportMinutesY, ...
                'supportWavelength', psfSupportWavelength, ...
                'zCoeffs', zCoeffs);
        end
        
    case 'Polans2015'
        % Polans optics
        for oiIndex = 1:oiNum
            %fprintf('Generating %s optics for eccentricity: %2.1f,%2.1f degs (um/deg):%2.1f\n', ...
            %    zernikeDataBase, oiSamplingGridDegs(oiIndex,1), oiSamplingGridDegs(oiIndex,2), obj.micronsPerDegree);
            targetEcc = oiSamplingGridDegs(oiIndex,:);
            
            [theOI, thePSF, psfSupportMinutesX, psfSupportMinutesY, psfSupportWavelength, zCoeffs] = ...
                PolansOptics.oiForSubjectAtEccentricity(subjectID, ...
                obj.whichEye, targetEcc, pupilDiamMM, obj.wave, obj.micronsPerDegree, ...
                'wavefrontSpatialSamples', wavefrontSpatialSamples, ...
                'subtractCentralRefraction', subtractCentralRefraction, ...
                'zeroCenterPSF', zeroCenterPSF, ...
                'flipPSFUpsideDown', flipPSFUpsideDown, ...
                'refractiveErrorDiopters', p.Results.refractiveErrorDiopters);
            
            oiEnsemble{oiIndex} = theOI;
            psfEnsemble{oiIndex} = struct(...
                'data', thePSF, ...
                'supportX', psfSupportMinutesX, ...
                'supportY', psfSupportMinutesY, ...
                'supportWavelength', psfSupportWavelength, ...
                'zCoeffs', zCoeffs);
        end
        
    end

end
