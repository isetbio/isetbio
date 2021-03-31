function [oiEnsemble, psfEnsemble] = oiEnsembleGenerate(obj, oiSamplingGridDegs, varargin)
    % Help
    if (ischar(oiSamplingGridDegs)) && (strcmp(oiSamplingGridDegs,'help'))
        doc('cMosaic.oiEnsembleGenerate');
        return;
    end
    
    % Parse input
    p = inputParser;
    p.addRequired('obj', @(x)(isa(x,'cMosaic')));
    p.addRequired('oiSamplingGridDegs', @(x)(isnumeric(x) && (size(x,2) == 2) && (all(isreal(x)))));
    p.addParameter('zernikeDataBase', 'Polans2015', @(x)(ismember(x, {'Polans2015'})));
    p.addParameter('subjectID', 6, @isscalar);
    p.addParameter('pupilDiameterMM', 3.0, @isscalar);
    p.parse(obj, oiSamplingGridDegs, varargin{:});

    oiSamplingGridDegs = p.Results.oiSamplingGridDegs;
    zernikeDataBase = p.Results.zernikeDataBase;
    pupilDiamMM = p.Results.pupilDiameterMM;
    subjectID = p.Results.subjectID;
    
    % Generate the oiEnsemble
    oiNum = size(oiSamplingGridDegs,1);
    oiEnsemble = cell(1, oiNum);
    psfEnsemble = cell(1, oiNum);
    
    switch (zernikeDataBase)
        case 'Polans2015'
            % Polans optics
            for oiIndex = 1:oiNum
                %fprintf('Generating %s optics for eccentricity: %2.1f,%2.1f degs (um/deg):%2.1f\n', ...
                %    zernikeDataBase, oiSamplingGridDegs(oiIndex,1), oiSamplingGridDegs(oiIndex,2), obj.micronsPerDegree);
                % determine average microns-per-deg for this eccentricity here
                [theOI, thePSF, psfSupportMinutesX, psfSupportMinutesY, psfSupportWavelength] = PolansOptics.oiForSubjectAtEccentricity(subjectID, ...
                    obj.whichEye, oiSamplingGridDegs(oiIndex,:), pupilDiamMM, obj.wave, obj.micronsPerDegree, ...
                    'wavefrontSpatialSamples', 301, ...
                    'subtractCentralRefraction', false);
                
                oiEnsemble{oiIndex} = theOI;
                psfEnsemble{oiIndex} = struct(...
                    'data', thePSF, ...
                    'supportX', psfSupportMinutesX, ...
                    'supportY', psfSupportMinutesY, ...
                    'supportWavelength', psfSupportWavelength);
            end
    end 
end
    