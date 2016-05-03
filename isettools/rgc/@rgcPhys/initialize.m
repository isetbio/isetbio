function obj = initialize(obj, ir, varargin)
% Initializes the rgcPhys object by loading a mosaic of GLM fits from an
% experiment in the Chichilnisky lab.
% 
% This function is only called by rgcPhys, which itself is only called by irPhys.
% 
%             rgcPhys = rgcPhys.initialize(rgc, varargin{:});   
% Inputs: 
%       rgc: an isetbio rgcPhys object
% 
% Outputs: the mosaic object, where each cell has a location, linear spatial
% and temporal receptive fields, a DC offest, a generator function, a
% post-spike filter, coupling filters if necessary, and empty fields for
% the linear, voltage and spiking responses.
% 
% 
% See also rgcPhys, irPhys.
% 
% (c) isetbio
% 09/2015 JRG%


%% Parse inputs
p = inputParser;
p.addRequired('obj');
p.addRequired('ir');
addParameter(p,'experimentID','2013-08-19-6',@ischar);
addParameter(p,'stimulusFit','WN',@ischar);
addParameter(p,'stimulusTest','NSEM',@ischar);
addParameter(p,'cellType','OnParasol',@ischar);


addParameter(p,'name','inner retina',@ischar);
addParameter(p,'species','unknown',@ischar);
% addParameter(p,'outersegment','os',@ischar);
addParameter(p,'eyeSide',    'left', @ischar);
addParameter(p,'eyeRadius',   4,     @isnumeric);
addParameter(p,'eyeAngle',    0,     @isnumeric);  % X-axis is 0, positive Y is 90

p.parse(obj,ir,varargin{:});

experimentID = p.Results.experimentID;
stimulusFit = p.Results.stimulusFit;
stimulusTest = p.Results.stimulusTest;
cellType =  p.Results.cellType;

obj.experimentID = experimentID;
obj.stimulusFit = stimulusFit;
obj.stimulusTest = stimulusTest;
obj.cellType = cellType;

%% Set defaults
obj.generatorFunction = @exp;
obj.numberTrials = 10;

% % Coupled experiment
% glmFitPath = '/Users/james/Documents/matlab/NSEM_data/';
% matFileNames = dir([glmFitPath '/ON*.mat']);

switch ieParamFormat(stimulusFit)
    case 'wn'        
        
        switch ieParamFormat(stimulusTest)
            case 'wn'
                glmFitPath = '/Users/james/Documents/matlab/akheitman/WN_mapPRJ/';
            case 'nsem'
                glmFitPath = '/Users/james/Documents/matlab/akheitman/WN_mapPRJ/Test_NSEM/';
        end
        
    otherwise % case 'NSEM'
        glmFitPath = '/Users/james/Documents/matlab/akheitman/NSEM_mapPRJ/';
end

switch ieParamFormat(cellType)
    case 'offparasol'
        matFileNames = dir([glmFitPath experimentID '/OFF*.mat']);        
    otherwise % case 'onparasol'
        matFileNames = dir([glmFitPath experimentID '/ON*.mat']);
end


% % % % % % 
% Loop through mat files and load parameters
for matFileInd = 1:2%118%length(matFileNames)
     
%     loadStr = sprintf('matFileNames(%d).name', matFileInd);
% %     eval(sprintf('load([glmFitPath %s])',loadStr))

%     fittedGLM = data.fittedGLM;

    cell = matFileNames(matFileInd).name(1:end-4);
    obj.cellID{matFileInd,1} = cell;
    load([glmFitPath experimentID '/' cell '.mat']);
    
%     
%     nameStr = eval(loadStr);
%     sind1 = strfind(nameStr,'_'); sind2 = strfind(nameStr,'.');
%     if isfield(fittedGLM.linearfilters,'Coupling')
% 
%         lookupIndex(matFileInd) = str2num(nameStr(sind1+1:sind2-1));
%     end
% %     lookupIndex(matFileInd) = 1205;
% %     fittedGLM = data.fittedGLM;
    
%     filterStimulus{matFileInd,1} = fittedGLM.linearfilters.Stimulus.Filter;
    obj.postSpikeFilter{matFileInd,1} = fittedGLM.linearfilters.PostSpike.Filter;
    if isfield(fittedGLM.linearfilters,'Coupling')

        obj.couplingFilter{matFileInd,1} = fittedGLM.linearfilters.Coupling.Filter;
    end
    
    obj.tonicDrive{matFileInd,1} = fittedGLM.linearfilters.TonicDrive.Filter;
    
    obj.sRFcenter{matFileInd,1} = fittedGLM.linearfilters.Stimulus.space_rk1;
    obj.sRFsurround{matFileInd,1} = 0*fittedGLM.linearfilters.Stimulus.space_rk1;
    obj.tCenter{matFileInd,1} = fittedGLM.linearfilters.Stimulus.time_rk1;
    obj.tSurround{matFileInd,1} = 0*fittedGLM.linearfilters.Stimulus.time_rk1;
    
    if isfield(fittedGLM.linearfilters,'Coupling')

        couplingMatrixTemp{matFileInd,1} = fittedGLM.cellinfo.pairs;
    end
    
    % NEED TO CHECK IF X AND Y ARE BEING SWITCHED INCORRECTLY HERE
    % figure; for i = 1:39; hold on; scatter(rgc2.mosaic{1}.cellLocation{i}(1), rgc2.mosaic{1}.cellLocation{i}(2)); end
    obj.cellLocation{matFileInd,1} = [fittedGLM.cellinfo.slave_centercoord.x_coord fittedGLM.cellinfo.slave_centercoord.y_coord];
    
%     % figure; imagesc(filterSpatial{matFileInd,1})
%     magnitude1STD = max(filterSpatial{matFileInd,1}(:))*exp(-1);
%     [cc,h] = contour(filterSpatial{matFileInd,1},[magnitude1STD magnitude1STD]);% close;
%     %         ccCell{rfctr} = cc(:,2:end);
%     cc(:,1) = [NaN; NaN];
%     spatialContours{matFileInd,1} = cc;
end

obj.rfDiameter = size(fittedGLM.linearfilters.Stimulus.Filter,1);
% if isfield(fittedGLM.linearfilters,'Coupling')
% 
% for matFileInd = 1:length(matFileNames)
% %     coupledCells = zeros(6,1);
%     for coupledInd = 1:length(couplingMatrixTemp{matFileInd,1})
%         coupledCells(coupledInd) = find(couplingMatrixTemp{matFileInd,1}(coupledInd)== lookupIndex);
%     end
%     obj.couplingMatrix{matFileInd,1} = coupledCells;    
%     
% end
% end
obj.couplingMatrix{1} = [17     3    11    34    12     9];

% % g = fittype('a*exp(-0.5*(x^2/Q1 + y^2/Q2)) + b*exp(-0.5*(x^2/Q1 + y^2/Q2))','independent',{'x','y'},'coeff',{'a','b','Q1','Q2'})
% 
%     ft = fittype( 'a*exp(-0.5*((x - x0)^2/Q1 + (y - y0)^2/Q2)) + b*exp(-0.5*((x - x0)^2/Q1 + (y - y0)^2/Q2)) + c0','independent',{'x','y'}, 'dependent', 'z', 'coeff',{'a','b','Q1','Q2','x0','y0','c0'});
%         opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     opts.StartPoint = [0.323369521886293 0.976303691832645];
%     
%     % Fit a curve between contrast level (x) and probability of correction
%     % detection.
%     % bootWeibullFit(stimLevels, nCorrect, nTrials, varargin)
%     % in computationaleyebrain/simulations/Pixel Visibility/ ...
%     [xsz,ysz] = size(srf1); 
%    [xc,yc] = meshgrid(1:xsz,1:ysz);
%     [fitresult, gof] = fit([xc(:),yc(:)], srf1(:), ft);, opts );
% 
% % Loop through mat files and plot contours
% figure; hold on;
% for matFileInd = 1:length(matFileNames)
%     plot(filterCenter{matFileInd,1}(1) + spatialContours{matFileInd,1}(1,2:end), filterCenter{matFileInd,1}(2) + spatialContours{matFileInd,1}(2,2:end))
%     
% end