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
addParameter(p,'fov',         8,     @isnumeric);
addParameter(p,'cellIndices',   0,     @isnumeric);
addParameter(p,'inputSize',[],           @isnumeric);
addParameter(p,'averageMosaic', 0,      @isnumeric);

p.parse(obj,ir,varargin{:});

experimentID = p.Results.experimentID;
stimulusFit = p.Results.stimulusFit;
stimulusTest = p.Results.stimulusTest;
cellType =  p.Results.cellType;
cellIndices = p.Results.cellIndices;
inputSize = p.Results.inputSize;

fov = p.Results.fov;
ecc = p.Results.eyeRadius;

averageFlag = p.Results.averageMosaic;

obj.experimentID = experimentID;
obj.stimulusFit = stimulusFit;
obj.stimulusTest = stimulusTest;
obj.cellType = cellType;

%% Set defaults
% obj.generatorFunction = @exp;
obj.numberTrials = 10;

% % Coupled experiment
% glmFitPath = '/Users/james/Documents/matlab/NSEM_data/';
% matFileNames = dir([glmFitPath '/ON*.mat']);

% switch ieParamFormat(stimulusFit)
%     case 'wn'        
%         
%         switch ieParamFormat(stimulusTest)
%             case 'wn'
%                 glmFitPath = '/Users/james/Documents/matlab/akheitman/WN_mapPRJ/';
%             case 'nsem'
%                 glmFitPath = '/Users/james/Documents/matlab/akheitman/WN_mapPRJ/Test_NSEM/';
%         end
%         
%     otherwise % case 'NSEM'
%         glmFitPath = '/Users/james/Documents/matlab/akheitman/NSEM_mapPRJ/';
% end



% RDT initialization
rdt = RdtClient('isetbio');

switch ieParamFormat(cellType)
    
    % RPE data set - need to put on RDT
    case 'onparasolrpe'
%         load('/Users/james/Documents/MATLAB/isetbio misc/rpeNora/mosaicGLM_RPE_onPar.mat')     

        rdt.crp('resources/data/rgc/rpe_dataset');
        data = rdt.readArtifact('mosaicGLM_RPE_onPar', 'type', 'mat');
        mosaicGLM = data.mosaicGLM;
    case 'offparasolrpe'
%         load('/Users/james/Documents/MATLAB/isetbio misc/rpeNora/mosaicGLM_RPE_offPar.mat')     

        rdt.crp('resources/data/rgc/rpe_dataset');
        data = rdt.readArtifact('mosaicGLM_RPE_offPar', 'type', 'mat');
        mosaicGLM = data.mosaicGLM;
    case 'onmidgetrpe'
%         load('/Users/james/Documents/MATLAB/isetbio misc/rpeNora/mosaicGLM_RPE_onMid.mat')     

        rdt.crp('resources/data/rgc/rpe_dataset');
        data = rdt.readArtifact('mosaicGLM_RPE_onMid', 'type', 'mat');
        mosaicGLM = data.mosaicGLM;
    case 'offmidgetrpe'
%         load('/Users/james/Documents/MATLAB/isetbio misc/rpeNora/mosaicGLM_RPE_offMid.mat')     

        rdt.crp('resources/data/rgc/rpe_dataset');
        data = rdt.readArtifact('mosaicGLM_RPE_offMid', 'type', 'mat');
        mosaicGLM = data.mosaicGLM;
 
    case {'onsbcrpe','sbcrpe'}
%         load('/Users/james/Documents/MATLAB/isetbio misc/rpeNora/mosaicGLM_RPE_onSBC.mat')     

        rdt.crp('resources/data/rgc/rpe_dataset');
        data = rdt.readArtifact('mosaicGLM_RPE_onSBC', 'type', 'mat');
        mosaicGLM = data.mosaicGLM;
    case 'offparasol'
%         matFileNames = dir([glmFitPath experimentID '/OFF*.mat']);        
%         load('/Users/james/Documents/MATLAB/isetbio misc/RDT uploads/mosaicGLM_WN_OFFParasol_2013_08_19_6.mat')

        rdt.crp('resources/data/rgc');
        data = rdt.readArtifact('mosaicGLM_WN_OFFParasol_2013_08_19_6', 'type', 'mat');
        mosaicGLM = data.mosaicGLM;
        
%         load('/Users/james/Documents/MATLAB/isetbio misc/RDT uploads/goodind_2013_08_19_6_OFFParasol.mat')

        rdt.crp('resources/data/rgc');                              
        data2 = rdt.readArtifact('goodind_2013_08_19_6_OFFParasol', 'type', 'mat');
        goodind = data2.goodind;
    otherwise % case 'onparasol'
%         matFileNames = dir([glmFitPath experimentID '/ON*.mat']);

        rdt.crp('resources/data/rgc');
        data = rdt.readArtifact('mosaicGLM_WN_ONParasol_2013_08_19_6', 'type', 'mat');
        mosaicGLM = data.mosaicGLM;        
        
        data2 = rdt.readArtifact('goodind_2013_08_19_6_ONParasol', 'type', 'mat');
        goodind = data2.goodind;
end

goodind = 1:length(mosaicGLM);


switch(averageFlag)
    
    % Load in a mosaic from a physiology experiment
    case 0
        
        if cellIndices ~= 0
            cellIndicesEval = cellIndices;
        else
            cellIndicesEval = [1:length(goodind)];
        end
        matFileCtr = 0;
        % % % % % %
        % Loop through mat files and load parameters
        % for matFileInd = 1:length(mosaicGLM)
        for matFileInd = cellIndicesEval%1:length(goodind)
            %     cell = matFileNames(matFileInd).name(1:end-4);
            %     obj.cellID{matFileInd,1} = cell;
            %     load([glmFitPath experimentID '/' cell '.mat']);
            matFileCtr = matFileCtr+1;
            obj.cellID{matFileCtr,1} = mosaicGLM{matFileInd}.cell_savename;
            %     obj.cellID{matFileInd,1} = mosaicGLM{goodind(matFileInd)}.cell_savename;
            
            if isfield(mosaicGLM{goodind(matFileInd)}.linearfilters,'PostSpike')
                obj.postSpikeFilter{matFileCtr,1} = mosaicGLM{goodind(matFileInd)}.linearfilters.PostSpike.Filter;
            else
                %         load('/Users/james/Documents/MATLAB/isetbio misc/rpeNora/psf1.mat')
                obj.postSpikeFilter{matFileCtr,1} = 0;%psf;
            end
            if isfield(mosaicGLM{goodind(matFileInd)}.linearfilters,'Coupling')
                
                obj.couplingFilter{matFileCtr,1} = mosaicGLM{goodind(matFileInd)}.linearfilters.Coupling.Filter;
            end
            
            
            
            switch ieParamFormat(cellType)
                case 'onparasolrpe'
                    obj.tonicDrive{matFileCtr,1} = 0;
                    obj.generatorFunction{matFileCtr,1} = mosaicGLM{goodind(matFileInd)}.model;
                otherwise
                    obj.tonicDrive{matFileCtr,1} = mosaicGLM{goodind(matFileInd)}.linearfilters.TonicDrive.Filter;
                    obj.generatorFunction{matFileCtr,1} = @exp;
            end
            
            obj.sRFcenter{matFileCtr,1} = mosaicGLM{goodind(matFileInd)}.linearfilters.Stimulus.space_rk1;
            obj.sRFsurround{matFileCtr,1} = 0*mosaicGLM{goodind(matFileInd)}.linearfilters.Stimulus.space_rk1;
            obj.tCenter{matFileCtr,1} = mosaicGLM{goodind(matFileInd)}.linearfilters.Stimulus.time_rk1;
            obj.tSurround{matFileCtr,1} = 0*mosaicGLM{goodind(matFileInd)}.linearfilters.Stimulus.time_rk1;
            
            if isfield(mosaicGLM{goodind(matFileInd)}.linearfilters,'Coupling')
                
                couplingMatrixTemp{matFileCtr,1} = mosaicGLM{goodind(matFileInd)}.cellinfo.pairs;
            end
            
            obj.cellLocation{matFileCtr,1} = [mosaicGLM{goodind(matFileInd)}.cellinfo.slave_centercoord.x_coord mosaicGLM{goodind(matFileInd)}.cellinfo.slave_centercoord.y_coord];
            
        end
        
        obj.rfDiameter = size(mosaicGLM{goodind(matFileInd)}.linearfilters.Stimulus.Filter,1);
        
        
    % Generate a mosaic that fully tiles the visual field 
    % with the average properties of an experimental mosaic,
    % shifted to a specified eccentricity.
    case 1
        
        
        switch ieParamFormat(cellType)
            
            % RPE data set - need to put on RDT
            
            case 'onparasolrpe'
                b = 10.7629; m = 18.9211;
                ecc0 = 10.9;
            case 'offparasolrpe'
                b = 68.3968; m = 0.85*18.9211;
                ecc0 = 10.9;
            case 'onmidgetrpe'
                b = -8.110; m = 10.7629;    
                ecc0 = 10.9;
            case 'offmidgetrpe'
                b = -8.110; m = 0.85*10.7629;   
                ecc0 = 10.9;
            case {'onsbcrpe','sbcrpe'}
                b = 70.2865; m = 15.8208;
                ecc0 = 10.9;
        end
        y0 = m*ecc0 + b;
        
        DFlimit = 4.5;
        y1 = max([DFlimit (m*ecc + b)]);
        
        yratio = y0/y1; 
        
        fov0 = 8;
        fovratio = fov/fov0;
        
        mosaicAverageGLM = mosaicAverage(mosaicGLM);
        
        cellSpacing = sqrt(mean(mosaicAverageGLM.sd.^2));
        
        stimRows = 40; stimCols = 80;
        
%         numberRows = floor(fovratio*yratio*stimRows/(2*cellSpacing));
%         numberCols = floor(fovratio*yratio*stimCols/(2*cellSpacing));
        
        numberCellsPerDegree0 = (stimCols/(2*cellSpacing))/fov0;
        numberCellsPerDegree = yratio*numberCellsPerDegree0;
        numberBipolarsPerDegree = inputSize(2)/fov;
        numberBipolarsPerCell = numberBipolarsPerDegree/numberCellsPerDegree;
        numberRows = floor(inputSize(1)/numberBipolarsPerCell);
        numberCols = floor(inputSize(2)/numberBipolarsPerCell);
        
%         numberPixelsPerCell0 = size(mosaicAverageGLM.linearfilters.Stimulus.space_rk1,1);
%         numberPixelsPerDegree0 = numberPixelsPerCell0 * numberCellsPerDegree0;
%         numberPixelsPerCell = numberPixelsPerDegree0 / numberCellsPerDegree;
        
        goodind = 1:numberRows*numberCols;
        
        matFileCtr = 0;
        for xi = 1:numberRows
            for yi = 1:numberCols
                matFileCtr = matFileCtr+1;
                
                obj.cellID{matFileCtr,1} = matFileCtr;
                
                if isfield(mosaicAverageGLM.linearfilters,'PostSpike')
                    obj.postSpikeFilter{matFileCtr,1} = mosaicAverageGLM.linearfilters.PostSpike.Filter;
                else
                    %         load('/Users/james/Documents/MATLAB/isetbio misc/rpeNora/psf1.mat')
                    obj.postSpikeFilter{matFileCtr,1} = 0;%psf;
                end
                if isfield(mosaicAverageGLM.linearfilters,'Coupling')
                    
                    obj.couplingFilter{matFileCtr,1} = mosaicAverageGLM.linearfilters.Coupling.Filter;
                end
                
                
                
                switch ieParamFormat(cellType)
                    case {'onparasolrpe','offparasolrpe','onmidgetrpe','offmidgetrpe','onsbcrpe','sbcrpe'}
                        obj.tonicDrive{matFileCtr,1} = 0;
                        obj.generatorFunction{matFileCtr,1} = mosaicAverageGLM.model;
                    otherwise
                        obj.tonicDrive{matFileCtr,1} = mosaicAverageGLM.linearfilters.TonicDrive.Filter;
                        obj.generatorFunction{matFileCtr,1} = @exp;
                end
                
                obj.sRFcenter{matFileCtr,1} = mosaicAverageGLM.linearfilters.Stimulus.space_rk1;
                % obj.sRFcenter{matFileCtr,1} = imresize(mosaicAverageGLM.linearfilters.Stimulus.space_rk1,[1+2*floor(numberBipolarsPerCell/2) 1+2*floor(numberBipolarsPerCell/2) ]);
                
                obj.sRFsurround{matFileCtr,1} = 0*mosaicAverageGLM.linearfilters.Stimulus.space_rk1;
                
                obj.tCenter{matFileCtr,1} = mosaicAverageGLM.linearfilters.Stimulus.time_rk1;
                obj.tSurround{matFileCtr,1} = 0*mosaicAverageGLM.linearfilters.Stimulus.time_rk1;
                
                if isfield(mosaicAverageGLM.linearfilters,'Coupling')
                    
                    couplingMatrixTemp{matFileCtr,1} = mosaicAverageGLM.cellinfo.pairs;
                end
                
                % obj.cellLocation{matFileCtr,1} = [(xi*2*cellSpacing + mod(yi,2)*cellSpacing) yi*2*cellSpacing];
                obj.cellLocation{matFileCtr,1} = [(xi*1*numberBipolarsPerCell + mod(yi,2)*0.5*numberBipolarsPerCell) yi*1*numberBipolarsPerCell];
                
            end
            
            obj.rfDiameter = size(mosaicAverageGLM.linearfilters.Stimulus.space_rk1,1);
        end
        
end