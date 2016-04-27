function obj = initialize(obj, rgc, cellTypeInd, varargin)
% intialize: a method of @rgcPhys that initializes the object
% following initialization by the superclass. This adds the generator
% function, the post spike filter and the coupling filters. This function
% is only called by rgcMosaicGLM, which itself is only called by rgcGLM or
% rgcCreate.
% 
%       rgcMosaicGLM.initialize(rgc, sensor, outersegment, varargin{:});
% 
% Inputs: 
%       rgc: an isetbio rgcGLM object
%       scene: an isetbio scene structure
%       sensor: an isetbio sensor structure
%       os: an isetbio outer segment structure
% 
% Outputs: the mosaic object with the generatorFunction, postSpikeFilter and
% couplingFilter properties set to appropriate values.
% 
% Example:
% 
%       obj.initialize(rgc, sensor, outersegment, varargin{:});
% 
% See also rgcCreate, rgcGLM, rgcMosaicGLM.
% 
% (c) isetbio
% 09/2015 JRG
 
%     % obj.generatorFunction = @erf;
%     obj.generatorFunction = @exp;
%     % obj.generatorFunction = @(x) 10*erf(x);
% 
%     obj.numberTrials = 10;
%     
%     obj.postSpikeFilter = buildPostSpikeFilter(.01);
%     
%     [obj.couplingFilter, obj.couplingMatrix] = buildCouplingFilters(obj, .01);


namesCellTypes = {'onParasol';'offParasol';'onMidget';'offMidget';'smallBistratified'};
obj.cellType = namesCellTypes{1};

obj.generatorFunction = @exp;

obj.numberTrials = 10;

glmFitPath = pwd;%'/Users/james/Documents/matlab/NSEM_data/';

% client = RdtClient('isetbio');
% % client.credentialsDialog();
% client.crp('resources/data/rgc');
% [data, artifact] = client.readArtifact('parasol_on_1205', 'type', 'mat');

glmFitPath = '/Users/james/Documents/matlab/NSEM_data/';
matFileNames = dir([glmFitPath '/ON*.mat']);

% Loop through mat files and load parameters
for matFileInd = 1:length(matFileNames)
     
    loadStr = sprintf('matFileNames(%d).name', matFileInd);
    eval(sprintf('load([glmFitPath %s])',loadStr))
    
    nameStr = eval(loadStr);
    sind1 = strfind(nameStr,'_'); sind2 = strfind(nameStr,'.');
    lookupIndex(matFileInd) = str2num(nameStr(sind1+1:sind2-1));
%     lookupIndex(matFileInd) = 1205;
%     fittedGLM = data.fittedGLM;
    
%     filterStimulus{matFileInd} = fittedGLM.linearfilters.Stimulus.Filter;
    obj.postSpikeFilter{matFileInd} = fittedGLM.linearfilters.PostSpike.Filter;
    obj.couplingFilter{matFileInd} = fittedGLM.linearfilters.Coupling.Filter;
    obj.tonicDrive{matFileInd} = fittedGLM.linearfilters.TonicDrive.Filter;
    
    obj.sRFcenter{matFileInd} = fittedGLM.linearfilters.Stimulus.space_rk1;
    obj.sRFsurround{matFileInd} = 0*fittedGLM.linearfilters.Stimulus.space_rk1;
    obj.tCenter{matFileInd} = fittedGLM.linearfilters.Stimulus.time_rk1;
    obj.tSurround{matFileInd} = 0*fittedGLM.linearfilters.Stimulus.time_rk1;
    
    couplingMatrixTemp{matFileInd} = fittedGLM.cellinfo.pairs;
    
    % NEED TO CHECK IF X AND Y ARE BEING SWITCHED INCORRECTLY HERE
    % figure; for i = 1:39; hold on; scatter(rgc2.mosaic{1}.cellLocation{i}(1), rgc2.mosaic{1}.cellLocation{i}(2)); end
    obj.cellLocation{matFileInd} = [fittedGLM.cellinfo.slave_centercoord.x_coord fittedGLM.cellinfo.slave_centercoord.y_coord];
    
%     % figure; imagesc(filterSpatial{matFileInd})
%     magnitude1STD = max(filterSpatial{matFileInd}(:))*exp(-1);
%     [cc,h] = contour(filterSpatial{matFileInd},[magnitude1STD magnitude1STD]);% close;
%     %         ccCell{rfctr} = cc(:,2:end);
%     cc(:,1) = [NaN; NaN];
%     spatialContours{matFileInd} = cc;
end

obj.rfDiameter = size(fittedGLM.linearfilters.Stimulus.Filter,1);

for matFileInd = 1:length(matFileNames)
    for coupledInd = 1:length(couplingMatrixTemp{matFileInd})
        coupledCells(coupledInd) = find(couplingMatrixTemp{matFileInd}(coupledInd)== lookupIndex);
    end
    obj.couplingMatrix{matFileInd} = coupledCells;    
    
end

% obj.couplingMatrix{1} = [17     3    11    34    12     9];

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
%     plot(filterCenter{matFileInd}(1) + spatialContours{matFileInd}(1,2:end), filterCenter{matFileInd}(2) + spatialContours{matFileInd}(2,2:end))
%     
% end