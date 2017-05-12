function ellipseParameters = ellipseGen(nRows, nCols, varargin)
% ellipseGen
% 
% Generate a number of [major axis, minor axis, rotation angle] ellipse
% parameter vectors for the spatial receptive fields of an RGC mosaic.
% 
% 3/2017 JRG (c) isetbio team

%%
p = inputParser;
p.addRequired('nRows'); % number rows of RGCs
p.addRequired('nCols'); % number cols of RGCs
p.addParameter('angleValues',[],@isnumeric);  % array of specific angle values
p.addParameter('axisValues',[],@isnumeric);   % array of specific axis values
p.addParameter('axisVariance',0.1,@isnumeric); % number for variance of major and minor axes
p.addParameter('ellipseParams',[],@isnumeric);% a single ellipseParams vector applied to every RGC 
p.parse(nRows,nCols,varargin{:});

nRows = p.Results.nRows;
nCols = p.Results.nCols;
angleValues = p.Results.angleValues;
axisValues = p.Results.axisValues;
axisVariance = p.Results.axisVariance;
ellipseParams =  p.Results.ellipseParams;

%% Build vectors

% Get values for major and minor axes if they are not already set
if isempty(axisValues)
    axisValues = (axisVariance*randn([nRows nCols 2]) + .5); zeroInd = axisValues<0; axisValues(zeroInd) = abs(axisValues(zeroInd));    
end

% Get angle values if they are not already set
if isempty(angleValues)
    angleValues = 180*(rand([nRows nCols])-.5);
end

% Build cell array for ellipse parameters out
ellipseParameters = cell(nRows,nCols);

% If the params are not specified by user with input, then generate them
if isempty(ellipseParams)
    
    for ii = 1:nRows
        for jj = 1:nCols
            ellipseParameters{ii,jj} = [squeeze(axisValues(ii,jj,:))',angleValues(ii,jj)];
        end
    end
    
else
    
    for ii = 1:nRows
        for jj = 1:nCols
            ellipseParameters{ii,jj} = [ellipseParams(1:2), ellipseParams(3)];
        end
    end
    
end

