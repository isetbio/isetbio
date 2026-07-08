classdef FakeConeMosaicForStrehlTest
    properties
        ExpectedZernikeDataBase
    end

    methods
        function obj = FakeConeMosaicForStrehlTest(expectedZernikeDataBase)
            obj.ExpectedZernikeDataBase = expectedZernikeDataBase;
        end

        function [oiEnsemble, psfEnsemble] = oiEnsembleGenerate(obj, oiPositionDegs, varargin)
            p = inputParser;
            p.addRequired('oiPositionDegs', @(x)(isnumeric(x) && isequal(size(x), [1 2])));
            p.addParameter('zernikeDataBase', '', @ischar);
            p.addParameter('subjectID', [], @isscalar);
            p.addParameter('pupilDiameterMM', [], @isscalar);
            p.addParameter('refractiveErrorDiopters', [], @isscalar);
            p.addParameter('noLCA', false, @islogical);
            p.addParameter('zeroCenterPSF', false, @islogical);
            p.addParameter('subtractCentralRefraction', false, @islogical);
            p.addParameter('wavefrontSpatialSamples', [], @isscalar);
            p.addParameter('upsampleFactor', [], @isscalar);
            p.addParameter('warningInsteadOfErrorForBadZernikeCoeffs', false, @islogical);
            p.parse(oiPositionDegs, varargin{:});

            assert(strcmp(p.Results.zernikeDataBase, obj.ExpectedZernikeDataBase), ...
                'FakeConeMosaicForStrehlTest:UnexpectedZernikeDataBase', ...
                'Unexpected zernikeDataBase value.');

            amplitude = localAmplitude(p.Results.subjectID, p.Results.refractiveErrorDiopters);
            oiEnsemble = {struct(...
                'oiPositionDegs', oiPositionDegs, ...
                'zernikeDataBase', p.Results.zernikeDataBase, ...
                'refractiveErrorDiopters', p.Results.refractiveErrorDiopters)};
            psfEnsemble = {struct(...
                'data', localPSFData(amplitude), ...
                'supportWavelength', [500 550 600])};
        end
    end
end

function amplitude = localAmplitude(subjectID, refractiveErrorDiopters)
if subjectID == 0
    amplitude = 2.0;
else
    amplitude = 1.0 - abs(refractiveErrorDiopters - 0.25) * 0.5;
end
end

function data = localPSFData(amplitude)
data = zeros(3, 3, 3);
data(2, 2, 1) = amplitude * 0.25;
data(2, 2, 2) = amplitude;
data(2, 2, 3) = amplitude * 0.5;
end
