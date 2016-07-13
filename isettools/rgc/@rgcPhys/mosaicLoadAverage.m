function obj = mosaicLoadAverage(obj, mosaicGLM, cellType, varargin)
% Load the parameters for an RGC mosaic measured in an experiment by the
% Chichilnisky Lab, find the average values of each paramter and store a
% mosaic with the average parameters in an isetbio object.
%
% JRG (c) 2016 isetbio team

%% Parse inputs
p = inputParser;
p.addRequired('obj');
p.addRequired('mosaicGLM');
p.addRequired('cellType');
addParameter(p,'cellIndices',   4,     @isnumeric);
addParameter(p,'goodind',    0,     @isnumeric);
addParameter(p,'params', [], @isstruct);
% addParameter(p,'ecc',   0.05,   @isnumeric);
% addParameter(p,'fov',   0.5,   @isnumeric);
p.parse(obj, mosaicGLM,cellType,varargin{:});

cellType = p.Results.cellType;
cellIndices = p.Results.cellIndices;
goodind  = p.Results.goodind;

params = p.Results.params;
% ecc = p.Results.ecc;
% fov = p.Results.fov;


experimentID = p.Results.params.experimentID;
stimulusFit  = p.Results.params.stimulusFit;
stimulusTest = p.Results.params.stimulusTest;
cellType     =  p.Results.params.cellType;
cellIndices  = p.Results.params.cellIndices;
inputSize    = p.Results.params.inputSize;
inputScale   = p.Results.params.inputScale;

fov = p.Results.params.fov;
ecc = p.Results.params.eyeRadius;
%% Calculate scaling parameters for mosaic at new eccentricity

% Scale the spatial RF size according to cell type
% Fit values from Dacey 2000 and Croner & Kaplan 1994
switch ieParamFormat(cellType)
    
    case {'onparasol','onparasolrpe','onparasolapricot'}
        b = 25;%68.3968/3;
        m = 18.9211;
        b2 = .0218; m2 = .0065;
        ecc0 = 10.9; % Do TEE conversion
    case {'offparasol','offparasolrpe','offparasolapricot'}
        b = 25;%68.3968/3;
        m = 0.85*18.9211;
        m = 0.85*18.9211;
        b2 = .0218; m2 = 0.85*.0065;
        ecc0 = 10.9;
    case {'onmidget','onmidgetrpe','onmidgetapricot'}
        b = -8.110; m = 10.7629;
        b2 = .0059; m2 = .0034;
        ecc0 = 10.9;
    case {'offmidget','offmidgetrpe','offmidgetapricot'}
        b = -8.110; m = 0.85*10.7629;
        b2 = .0059; m2 = 0.85*.0034;
        ecc0 = 10.9;
    case {'sbc','onsbcrpe','sbcrpe','sbcapricot'}
        b = 35;%70.2865;
        m = 20;%15.8208;
        b2 = 70.2865; m2 = 15.8208;
        ecc0 = 10.9;
end

% Calculate conversion ratio
y0 = m*ecc0 + b;
y02 = m2*ecc0 + b2;

% Set limit for nearing zero eccentricity
DFlimit = 1.1*4.5;
y1 = max([DFlimit (m*ecc + b)]);
y12 = max([0 (m2*ecc + b2)]);

% Scale RF size
yratio = y0/y1;
y2ratio = y02/y12;

% Scale field of view
fov0 = 8;
fovratio = fov/fov0;

% Compute a mosaic of average parameters
mosaicAverageGLM = mosaicAverage(mosaicGLM);

% Get average 1 sd spacing
cellSpacing = sqrt(mean(mosaicAverageGLM.sd.^2));

% Set parameters from physiology experiments
stimRows = 40; stimCols = 80;

% Build a hexagonal lattice in periphery and square in fovea
hexFlag = 1;

% Calculate cell density on new mosaic
numberCellsPerDegree0 = (stimCols/(2*cellSpacing))/fov0;
numberCellsPerDegree = yratio*numberCellsPerDegree0;

% InputScale incorporates number of cones per bipolar
numberBipolarsPerDegree = inputScale*inputSize(2)/(fov);
numberBipolarsPerCell = (numberBipolarsPerDegree/numberCellsPerDegree);

% Set one cone to one bipolar in high density case
if numberBipolarsPerCell < 2
    numberBipolarsPerCell = 1;
    hexFlag = 0;
end

% Calculate number of cells in new mosaic
numberRows = floor(1*inputSize(1)/numberBipolarsPerCell);
numberCols = floor(1*inputSize(2)/numberBipolarsPerCell);

% Get original cell spacing and compute new spacing
cellSpacingRF = size(mosaicAverageGLM.linearfilters.Stimulus.space_rk1,1);
numberCellsPerDegreeRF0 = (stimCols/(2*cellSpacingRF))/fov0;
numberCellsPerDegreeRF = yratio*numberCellsPerDegreeRF0;

numberBipolarsPerDegreeRF = inputScale*inputSize(2)/(fov);
numberBipolarsPerCellRF = numberBipolarsPerDegreeRF/numberCellsPerDegreeRF;

goodind = 1:numberRows*numberCols;

matFileCtr = 0;
for xi = 1:numberRows
    for yi = 1:numberCols
        matFileCtr = matFileCtr+1;
        
        obj.cellID{matFileCtr,1} = matFileCtr;
        
        % Post spike filter
        if isfield(mosaicAverageGLM.linearfilters,'PostSpike')
            obj.postSpikeFilter{matFileCtr,1} = mosaicAverageGLM.linearfilters.PostSpike.Filter;            
        else
            obj.postSpikeFilter{matFileCtr,1} = 0;
        end
        
        % Coupling filters
        if isfield(mosaicAverageGLM.linearfilters,'Coupling')
            obj.couplingFilter{matFileCtr,1} = mosaicAverageGLM.linearfilters.Coupling.Filter;
            couplingMatrixTemp{matFileCtr,1} = mosaicAverageGLM.cellinfo.pairs;
        end
        
        % 1 SD magnitude at which to draw contours
        obj.rfDiaMagnitude = inputScale*numberBipolarsPerCell;
        
        % Tonic drive and generator function
        switch ieParamFormat(cellType)
            case {'onparasolrpe','offparasolrpe','onmidgetrpe','offmidgetrpe','onsbcrpe','sbcrpe'}
                obj.tonicDrive{matFileCtr,1} = 0;
                obj.generatorFunction{matFileCtr,1} = mosaicAverageGLM.model;
            case{'onparasolapricot','offparasolapricot','onmidgetapricot','offmidgetapricot','onsbcapricot','sbcapricot'}
                obj.tonicDrive{matFileCtr,1} = 0;
                obj.generatorFunction{matFileCtr,1} = @exp;
            otherwise
                obj.tonicDrive{matFileCtr,1} = mosaicAverageGLM.linearfilters.Stimulus.tonicDrive;
                obj.generatorFunction{matFileCtr,1} = @exp;
        end
        
        % Resize spatial RF
        rf1 = imresize(mosaicAverageGLM.linearfilters.Stimulus.space_rk1,[1+2*floor(numberBipolarsPerCellRF/2) 1+2*floor(numberBipolarsPerCellRF/2) ]);
        rf2 = mosaicAverageGLM.linearfilters.Stimulus.space_rk1;
        rf1rs = (rf1.*sum(rf2(:))/sum(rf1(:)));
        obj.sRFcenter{matFileCtr,1} = rf1rs;
      
        % For rk1 fits, set surround RF to zero
        obj.sRFsurround{matFileCtr,1} = zeros(size(rf1rs));
        
        % Calculate average temporal response
        obj.tCenter{matFileCtr,1} = mosaicAverageGLM.linearfilters.Stimulus.time_rk1;
        obj.tSurround{matFileCtr,1} = 0*mosaicAverageGLM.linearfilters.Stimulus.time_rk1;        
        
        % Build mosaic locations
        obj.cellLocation{matFileCtr,1} = (1/1).*[(xi*1*numberBipolarsPerCell + hexFlag*mod(yi,2)*0.5*numberBipolarsPerCell) yi*1*numberBipolarsPerCell];
        
    end
    
    % Set spatial RF diameter for all cells
    obj.rfDiameter = size(obj.sRFcenter{matFileCtr,1},1);
end