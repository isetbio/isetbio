classdef rgcLayer < handle
    %Defines the rgc layer class
    %
    %   rgcL = rgcLayer;
    %
    % This class manages all the layer-specific parameters from an rgc
    % parameters class.
    %
    % When creating a layer, you can create several alternative defaults, (see
    % rgcParameters)
    %
    % To see a list of the properties use
    %
    %    properties('rgcLayer')   or
    %    methodsview('rgcLayer')
    %
    % I wonder if layers should be a subclass of the rgcParameters class? (BW)
    %
    % To read the parameter from a class we invoke
    %   layer.get('Property')
    %
    % The public properties are
    % (c) 2010 Stanford Synapse Team
    
    % TODO:
    %  Deal with the few functions all the way at the bottom.
    
    %% Public Properties
    properties (GetAccess = 'public', SetAccess = 'public')
        
        %% Book keeping
        name = '';
        
        % RGC parent simulation.  Explain how this gets set.
        parent;
        
        %% Electrical properties
        vSwing = 0.1;           % 100 millivolts is default
        rgcvTimeSeries = [];    % rgc voltages; nRGC x n time points
        rgcVoltThresh  = 0.08;  % The voltage required for a spike.  Must be less than vSwing
        
        %% Time series properties
        currentSpkTS  = [];   % The spike time series  (uint8s at Time samples)
        currentLinTS  = [];   % linear TS (cone-driven); nRGC x nTimePoints?
        RFcomponentTS = [];         % cone driven inputs before combining center and surround;
        % cell array, 1 x nRGC_Components;
        % each cell array is nRGC x nTimePoints
        
        
        %% Topology parameters
        currentConnec = [];   % connection matrix for
        % I changed the distanceFunction to
        % wRange*exp(-D)+rand(size(D))*nRange;
        wRange = 1;
        nRange = 0.02;
        cutoff = .05;   % See explanations below
        
        %% Set spatial RF properties
        cellSpacing = 1.5;         % dS between cells in um. 
        layerCenter = [0 0];       % center of the layer in um
        RFcov       = {[5 0; 0 5] [10 0; 0 10]};    % Spatial covariance matrix in um^2
        
        % [sxx sxy; syx syy] can be also set using:
        % rgcP.set('RF Covariance Matrix',{[sxx syy theta] [sxx syy theta]}, layerNumber);
        %RFsCoeffs   = [1; -.8];    % Center and surround weights
        
        % Let's make them balanced just to see what happens
        RFsCoeffs   = [1; -1];
        
        % How do the different cone classes contribute to the layer?
        coneWeights = ones(2,3); 
               
        % Create response parameters
        % These should all go away along with the corresponding buildXXX
        % functions below.
        %         centerBDiv   = [32 16];     % center.b = trDur ./ centerBDiv;
        %         centerF      = [1 .5];
        %         centerNormF  = .001;
        %         centerGamma  = [5 5];
        %
        %         surroundBDiv   = [28 12];   %surround.b = trDur ./ surroundBDiv;
        %         surroundF      = [1 .5];
        %         surroundNormF  = .001;
        %         surroundGamma  = [5 5];
        
        % Temporal response of the center and surround
        centerTR   = [];
        surroundTR = [];

        % The coupling function (i.e. the feedback function of one neuron
        % on its neighbors) fbTR and cpTR are trigerred only when a neuron
        % spikes
        % Feedback temporal response: corresponds to the feedback
        % function of one neuron on itself
        %  Set in rgcLayer call by 
        %  cpTR = layerInit('fbtr');
        fbTR = [];
        hasFeedback = 1;

        % Coupling temporal response
        %  Set in rgcLayer call by 
        %  cpTR = layerInit('cptr');
        cpTR = [];
        hasCoupling = 1;
        
    end
    
    %% No function - Not dependent because default is false.
    properties (GetAccess = 'public', SetAccess = 'private')
        overrideSizeDefault = 0;    % must be 0 or 1, is it overriden?
        overridenSize;              % in case the total size is overriden
        
        % measured in um.
        defaultCutoff = 0.05; % this is the initial cutoff value.
        
    end
    
    %% Properties defined from other parameters (dependent)
    %  These require a get function.
    properties(Dependent = true, SetAccess = 'private')
        % Basic properties
        gridSize;    % size of the layer measured in cells
        cellGrid     % Grid giving the locations of the cells (um)
        cellLoc      % Cell Locations as (X(:),Y(:)) positions in (um)
        cellConePos  % Cell positions with respect to the cone indices
        
        dT;         % Time step
        
        RF;         % RGC Receptive field
        
        % regarding center and surround, see twoGammaResp.m to have more
        % detail on their effect. These are used to compute inTR
        inTR;       % input temporal response
        % inTR has 2 dimensions as
        % inTR(:,1) represents the input temporal corresponding to the
        % center component of the RF, and inTR(:,2) corresponds to the
        % surround component of the RF.
        
        % Function used to compute the feedback temporal response
        
        trDur;      % Input response function duration (ms)
        
    end
    
    
    %% Properties calculated from defined parameters (dependent)
    methods
        function gs = get.gridSize(obj)
            % Number of cells (row,col)
            % The number of RGCs is determined by the number of cones.  If
            % you want to experiment with a small grid size, use one of the
            % small scenes that sweeps out a small field of view.
            units = 'um';
            % The obj is the layer.
            
            % row/col dimensions for the cells in the grid.
            % The spatial coverage of the entire cone array (um)
            rgcP = obj.parent;
            
            % Total size of the cone image in um
            coneImageSize = rgcP.get('cone image size',units);
            %cs = rgcP.get('cone grid size')*rgcP.get('coneSpacing');
            
            % We calculate how many cells are needed to cover the full cone
            % array. cellSpacing is in microns.  Should get a get with
            % units.
            % Cell center spacing in um
            gss = coneImageSize / layerGet(obj,'cell Spacing',units);
            
            gs = ceil(gss);
            
        end
        
        function cg = get.cellGrid(obj)
            % The grid locations in um
            gs = obj.gridSize;
            newCenter = (obj.overrideSizeDefault)*obj.layerCenter;
            cg{1} = (1:gs(1))*obj.cellSpacing;
            cg{1} = cg{1}-mean(cg{1})+newCenter(1);
            cg{2} = (1:gs(2))*obj.cellSpacing;
            cg{2} = cg{2}-mean(cg{2})+newCenter(2);
        end
        
        function cl = get.cellLoc(obj)
            % The grid locations as a meshgrid (in um)
            cg = obj.cellGrid;
            [X Y] = meshgrid(cg{1}', cg{2});
            cl = [X(:) Y(:)];
        end
        
        function pos = get.cellConePos(obj)
            % The position of the cells with respect to the cone indices
            rgcP = obj.parent;
            s = round(obj.get('cell spacing') / rgcP.get('cone spacing'));
            pos = cell(2,1);
            sz = rgcP.get('cone grid size');
            pos{1} = 1:s:sz(1);
            pos{2} = 1:s:sz(2);
        end
        
        % This is the time step in ms
        function dT = get.dT(obj)
            dT = obj.parent.get('dT');
        end
        
        function td = get.trDur(obj)
            td = obj.parent.get('trDur');
        end
        
        % Temporal response for the center
        function impulseResp = get.centerTR(obj)
            % center = buildCenter(obj);
            impulseResp = obj.centerTR;
        end
        
        % Temporal response for the surround
        function impulseResp = get.surroundTR(obj)
            % surround = buildSurround(obj);
            impulseResp = obj.surroundTR;
        end
        
        % This is a temporal response to the input.
        % Two  responses, one for the center and another the surround.
        function inTR = get.inTR(obj)
            % The twoGammaResp doesn't return vectors that sum to one. So,
            % we do that here.  Maybe we should force twoGammaResp to
            % return IRFs that sum to 1.
            % rgcP = obj.get('parent');
                        
            % The input (stimulus-driven)
            centerIRF   = obj.get('center tr');
            surroundIRF = obj.get('surround tr');
            
            inTR = [centerIRF(:) surroundIRF(:)];
            % plot(inTR)
        end
                
        function rf = get.RF(obj)
            % The returns a matrix in which rf(:,:,ii) is the RF of the ith
            % component for this layer.  If you want the sum, you should
            % use layer.get('rf sum')
            rf = rgcComputeRFs(obj);
            
        end
    end
    
    %% Class constructor
    methods
        function obj = rgcLayer(papa)
            if notDefined('papa')
                warning('This layer does not depend on any network you really should not do that.');
            else
                %papa should be in the class rgcParameters
                if (~isequal(class(papa),'rgcParameters'))
                    error('invalid kind of parent');
                else
                    obj.parent = papa;
                end
            end
            
            % Center and surround temporal responses 
            % These are now stored as time series (dt steps).
            obj.centerTR   = layerInit('center tr',papa);
            obj.surroundTR = layerInit('surround tr',papa);
            
            % Feedback and coupling
            obj.fbTR = layerInit('fb tr',papa);
            obj.cpTR = layerInit('cp tr',papa);
            
        end
    end
    
    %% basics (get/set, Display, ...) methods
    methods
        function val = get(obj,fieldName)
            val = layerGet(obj,fieldName);
        end
        
        function obj = set(obj,fieldName,value)
            obj = layerSet(obj,fieldName,value);
        end
        
        function Display(obj,fieldName)
            layerSet(obj,fieldName);
        end
        
        function resetCurrentState(obj)
            obj.rgcvTimeSeries   = [];
            obj.currentSpkTS  = [];
            obj.currentLinTS  = [];
            obj.currentConnec = [];
            % obj.currentSpikeIntegral = [];
        end
        
        function newObj = Copy(obj)
            newObj = rgcLayer(obj.parent);
            % public fields
            rgcCopy(obj,newObj);
            % private fields
            newObj.overrideSizeDefault = obj.overrideSizeDefault;
            newObj.overridenSize = obj.overridenSize;
        end
        
        function bord = getBorderSize(obj) % in um
            rf = obj.RF;
            rfS = size(rf);
            cS = obj.parent.coneSpacing;
            bord = [rfS(1) rfS(2)]*cS/2; % only half of the RF on each side of the cell
        end
        
        function resetSize(obj)
            obj.overrideSizeDefault = 0;
            obj.overridenSize = [];
            obj.layerCenter = [0 0];
        end
        
        function set.RFcov(obj,value)
            % Must have a covariance matrix for the center and the
            % surround.
            if ~iscell(value) || length(value) ~= 2
                error('RFcov must be a 1x2 cell');
            end
            for ii = 1:2
                if ~isCovarianceMatrix(value{ii})
                    error('RFcov{%d} must be 2x2, symmetric positive definite',ii)
                end
            end
            obj.RFcov = value;
        end
        
        
    end
end

%% Other functions - Why are they here? Maybe they should be up above?

% Initialization of layer variables
function val = layerInit(param,rgcP,varargin)
% Initialize temporal response parameters for feedback and coupling

% Default parameter values used at initialization
tShift = .05;
linF   = .02;

param = stringFormat(param);
switch(param)
    case {'fbtr','feedbacktemporalresponse'}
        % Feedback temporal response
        t = rgcP.get('tr samples','ms');
        
        % cf. with Pillow et al. JNS 2005, Figure 3ab, right panel
        fbDip = zeros(size(t));  fbDip(1:2) = -1*linF;
        fbSlow = zeros(size(t)); fbSlow(3:6) = linF; fbSlow(7:12) = -linF/2;
        % Could blur the fbSlow and fbDip here
        val = fbDip + fbSlow;
        % vcNewGraphWin; plot(val)
                
    case {'cptr','couplingtemporalresponse'}
        % Coupling temporal response
        t = rgcP.get('tr samples','ms');
        val   = (1 ./ (t .* tShift)) * linF;
        
    case 'centertr'  % Center impulse response function
        t = rgcP.get('tr samples','sec');
        
        % It would be best to specify the parameters in terms of first
        % peak, and total duration and size of negative compared to
        % positive. For now, we have this.
        tau = 0.015; n = 2;
        val   = gammaPDF(t,tau,n);
        
    case 'surroundtr'% Surround impulse response function
        t = rgcP.get('tr samples','sec');
        
        % The surround is smaller than the center
        tau = 0.025; n = 2;
        val   = 0.5*gammaPDF(t,tau,n);

    otherwise
        error('Unknown parameters %s\n',param);
end

end

% Should buildCenter and the others return the whole spatial array?
% function center = buildCenter(obj)
% % Temporal impulse response of center
% center = struct();
% center.b      = obj.trDur ./ obj.centerBDiv;
% center.f      = obj.centerF;
% center.normF  = obj.centerNormF;
% center.gamma  = obj.centerGamma;
% end

% Default parameters
% function surround = buildSurround(obj)
% % Temporal impulse response of surround
% surround = struct();
% surround.b      = obj.trDur ./ obj.surroundBDiv;
% surround.f      = obj.surroundF;
% surround.normF  = obj.surroundNormF;
% surround.gamma  = obj.surroundGamma;
% end

% Test if a matrix is positive definite covariance (2x2)
function res = isCovarianceMatrix(R)
res = 1;
res = res && all(size(R) == [2 2]);
res = res && (abs(R(1,2)-R(2,1)) < 1e-6);
res = res && (det(R) > 0);
res = res && (trace(R) > 0);
end
