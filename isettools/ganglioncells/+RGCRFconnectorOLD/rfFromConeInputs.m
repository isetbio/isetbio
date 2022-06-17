function [spatialSupportXY, RF, theInputConeLineSpreadFunctionsXY] = rfFromConeInputs(...
            theConePositions, theConeSpacings, theConeWeights, varargin)
% Generate the 2D RF of an RGC from its input cones
%
% Syntax:
%   [spatialSupportXY, RF, ...
%    theInputConeLineSpreadFunctionsXY] = ...
%         RGCRFconnector.rfFromConeInputs(...
%            theConePositions, theConeSpacings, theConeWeights, varargin)
%
% Description:
%   Generate the 2D RF of an RGC from its input cones
%
% Inputs:
%    theConePositions   - Vector with indices of the input cones
%    theConeSpacings    - Vector with the spacing of the input cones
%    theConeWeights     - Vector with the pooling weights of the input cones
%
%
% Outputs:
%    spatialSupportXY                   - [N x 2] matrix of RF spatial support (x/y)
%    RF                                 - [N x N] matrix of RF
%    theInputConeLineSpreadFunctionsXY  - {2 x inputsNum] matrix of the line spread
%                                         functions of the input cone apertures along 
%                                         the X- and the Y-axes
%
% Optional key/value pairs
%   'xSupport'                  - Custom spatial support (X)
%   'ySupport'                  - Custom spatial support (Y)
%   'spatialSupportSamples'     - Custom spatial support samples
%
% History:
%   5/11/2022       NPC     Wrote it
%

    p = inputParser;
    p.addOptional('xSupport', [], @(x)(isempty(x) || (isnumeric(x))));
    p.addOptional('ySupport', [], @(x)(isempty(x) || (isnumeric(x))));
    p.addOptional('spatialSupportSamples', 60, @(x)(isempty(x) || (isnumeric(x))));
    p.parse(varargin{:});
    xSupport = p.Results.xSupport;
    ySupport = p.Results.ySupport;
    spatialSupportSamples = p.Results.spatialSupportSamples;

    % Compute spatial support
    xSep = max(theConeSpacings);

    if (isempty(xSupport))
        xx = theConePositions(:,1);
        xSupport = linspace(min(xx)-xSep,max(xx)+xSep,spatialSupportSamples);
    end

    if (isempty(ySupport))
        yy = theConePositions(:,2);
        ySupport = linspace(min(yy)-xSep,max(yy)+xSep,spatialSupportSamples);
    end
    [X,Y] = meshgrid(xSupport, ySupport);
    spatialSupportXY(:,1) = xSupport(:);
    spatialSupportXY(:,2) = ySupport(:);
    RF = zeros(size(X));

    % Output
    theConeInputLineSpreadFunctionXY = cell(2, numel(theConeWeights));

    for iCone = 1:numel(theConeWeights)
        % Characteristic radius of the input cone
        rC = 0.204*sqrt(2.0)*theConeSpacings(iCone);
        % Compute aperture2D x weight
        XX = X-theConePositions(iCone,1);
        YY = Y-theConePositions(iCone,2);
        theConeAperture2D = theConeWeights(iCone) * exp(-(XX/rC).^2) .* exp(-(YY/rC).^2);
        % Accumulate 2D apertures
        RF = RF + theConeAperture2D;
        % 1D Line spread functions (XY) for each cone aperture
        theInputConeLineSpreadFunctionsXY{1,iCone} = sum(theConeAperture2D,1);
        theInputConeLineSpreadFunctionsXY{2,iCone} = sum(theConeAperture2D,2);
    end
    
end