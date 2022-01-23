function ellipseParameters = ellipseGen(nRows, nCols, varargin)
% Generate a number of ellipse parameter vectors for a RGC mosaic.
%
% Syntax:
%   ellipseParameters = ellipseGen(nRows, nCols, [varargin])
%
% Description:
%    Generate a number of [major axis, minor axis, rotation angle] ellipse
%    parameter vectors for the spatial receptive fields of an RGC mosaic.
%
% Inputs:
%    nRows             - Numeric. The number of rows of RGCs
%    nCols             - Numeric. The number of columns of RGCs
%
% Outputs:
%    ellipseParameters - Cell. A cell array of ellipse parameters for the
%                        retinal ganglion cells.
%
% Optional key/value pairs:
%    angleValues       - Array. An array of specific angle values. Default
%                        value is [].
%    axisValues        - Array. An array of specific axis values. Default
%                        value is [].
%    axisVariance      - Numeric. The number for the variance of major and
%                        minor axes. Default value is 0.1.
%    ellipseParams     - Array. A single ellipse parameters vector applied
%                        to every RGC. Default value is [].
%

% History:
%    03/XX/17  JRG  (c) isetbio team
%    06/04/19  JNM  Documentation pass

%%
p = inputParser;
p.KeepUnmatched = true;
p.addRequired('nRows'); % number rows of RGCs
p.addRequired('nCols'); % number cols of RGCs
p.addParameter('angleValues', [], @isnumeric);    % Specific angles array
p.addParameter('axisValues', [], @isnumeric);     % Specific axis array
p.addParameter('axisVariance', 0.1, @isnumeric);  % Major/Minor axis var.
% Single ellipseParams vector applied to every RGC 
p.addParameter('ellipseParams', [], @isnumeric);
p.parse(nRows, nCols, varargin{:});

nRows = p.Results.nRows;
nCols = p.Results.nCols;
angleValues = p.Results.angleValues;
axisValues = p.Results.axisValues;
axisVariance = p.Results.axisVariance;
ellipseParams =  p.Results.ellipseParams;

%% Build vectors
% Get values for major and minor axes if they are not already set
if isempty(axisValues)
    axisValues = (axisVariance * randn([nRows nCols 2]) + 1);
    zeroInd = axisValues < 0;
    axisValues(zeroInd) = abs(axisValues(zeroInd));    
end

% Get angle values if they are not already set
if isempty(angleValues)
    angleValues = 180 * (rand([nRows nCols]) - .5);
end

% Build cell array for ellipse parameters out
ellipseParameters = cell(nRows, nCols);

% If the params are not specified by user with input, then generate them
if isempty(ellipseParams)
    for ii = 1:nRows
        for jj = 1:nCols
            axisValues(ii, jj, :) = axisValues(ii, jj, :) ./ ...
                norm(squeeze([axisValues(ii, jj, :)]));
            ellipseParameters{ii, jj} = ...
                [squeeze(axisValues(ii, jj, :))', angleValues(ii, jj)];
        end
    end
else
    ellipseParams(1:2) = ellipseParams(1:2) ./ ...
        norm(squeeze([ellipseParams(1:2)]));
    for ii = 1:nRows
        for jj = 1:nCols
            ellipseParameters{ii, jj} = ...
                [ellipseParams(1:2), ellipseParams(3)];
        end
    end
end
