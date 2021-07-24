function [oiEnsemble, psfEnsemble, zCoeffs] = oiEnsembleGenerate(obj, oiSamplingGridDegs, varargin)
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
    p.addParameter('subjectID', 6, @isscalar);
    p.addParameter('pupilDiameterMM', 3.0, @isscalar);
    p.addParameter('wavefrontSpatialSamples', 301, @isscalar);
    p.addParameter('subtractCentralRefraction', false, @islogical);
    p.addParameter('zeroCenterPSF', true, @islogical);
    p.addParameter('flipPSFUpsideDown', true, @islogical);
    p.parse(obj, oiSamplingGridDegs, varargin{:});

    oiSamplingGridDegs = p.Results.oiSamplingGridDegs;
    zernikeDataBase = p.Results.zernikeDataBase;
    pupilDiamMM = p.Results.pupilDiameterMM;
    subjectID = p.Results.subjectID;
    subtractCentralRefraction = p.Results.subtractCentralRefraction;
    wavefrontSpatialSamples = p.Results.wavefrontSpatialSamples;
    zeroCenterPSF = p.Results.zeroCenterPSF;
    flipPSFUpsideDown = p.Results.flipPSFUpsideDown;
    
    % Generate the oiEnsemble
    oiNum = size(oiSamplingGridDegs,1);
    oiEnsemble = cell(1, oiNum);
    psfEnsemble = cell(1, oiNum);
    
    switch (zernikeDataBase)
        
        case 'Artal2012'
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
                    'flipPSFUpsideDown', flipPSFUpsideDown);
                
                if (isempty(theOI))
                    error('Bad Zernike coefficents for the %s of Artal subject %d. Choose another subject/eye', obj.whichEye, subjectID);
                end
                
                oiEnsemble{oiIndex} = theOI;
                psfEnsemble{oiIndex} = struct(...
                    'data', thePSF, ...
                    'supportX', psfSupportMinutesX, ...
                    'supportY', psfSupportMinutesY, ...
                    'supportWavelength', psfSupportWavelength);
            end
            
        case 'Polans2015'
            % Polans optics
            for oiIndex = 1:oiNum
                %fprintf('Generating %s optics for eccentricity: %2.1f,%2.1f degs (um/deg):%2.1f\n', ...
                %    zernikeDataBase, oiSamplingGridDegs(oiIndex,1), oiSamplingGridDegs(oiIndex,2), obj.micronsPerDegree);
                targetEcc = oiSamplingGridDegs(oiIndex,:);
                
                [theOI, thePSF, psfSupportMinutesX, psfSupportMinutesY, psfSupportWavelength, zCoeffs] = PolansOptics.oiForSubjectAtEccentricity(subjectID, ...
                    obj.whichEye, targetEcc, pupilDiamMM, obj.wave, obj.micronsPerDegree, ...
                    'wavefrontSpatialSamples', wavefrontSpatialSamples, ...
                    'subtractCentralRefraction', subtractCentralRefraction, ...
                    'zeroCenterPSF', zeroCenterPSF, ...
                    'flipPSFUpsideDown', flipPSFUpsideDown);
                
                oiEnsemble{oiIndex} = theOI;
                psfEnsemble{oiIndex} = struct(...
                    'data', thePSF, ...
                    'supportX', psfSupportMinutesX, ...
                    'supportY', psfSupportMinutesY, ...
                    'supportWavelength', psfSupportWavelength);
            end
    end 
end
    