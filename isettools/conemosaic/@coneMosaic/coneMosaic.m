classdef coneMosaic < hiddenHandle
    % Create a cone mosaic class
    %
    %   cMosaic =  coneMosaic('pigment', photoPigment, 'os', os);
    %
    % The cone mosaic defines an array of cones. The individual cones have
    % absorption properties defined by cMosaic.pigment. The computation
    % from absorptions to photocurrent is defined by cMosaic.os
    %
    % HJ/JRG/BW ISETBIO Team, 2016
    
    properties  % public properties
        name                % the name of the object
        
        pigment;            % Cone class object, contain single cone property
        macular;            % Macular class object
        os;                 % Outersegment properties
        
        center;             % (x,y) center position of patch in meters
        
        pattern;            % Pattern of K-LMS cones in the mosaic
        patternSampleSize;  % Separation between K-LMS pattern samples; for rectangular grid mosaics, 
                            % this is set to the pigment width/height, i.e., the actual cone separation;
                            % For hexagonal grid mosaics (instances of the coneMosaicHex class), 
                            % this becomes the separation between the rect grid nodes over which the 
                            % lower resolution hex grid is sampled (NC)
        
        integrationTime;    % Cone temporal integration time in secs (50 ms default)
        emPositions;        % Eye movement positions in number of cones.
                            % The length of this property controls number of
                            % frames to be computed
        noiseFlag;          % To control which noise is included
        hdl;                % handle of the gui window
    end
    
    %     properties (SetObservable, AbortSet)
    %         sampleTime;         % Time step for em and os computation. In the
    %                             % os this is called timeStep. Default is 1 ms
    %     end
    
    properties (GetAccess=public, SetAccess=public) % public temporarilly
        absorptions;    % The spatial array of cone absorptions
    end
    
    properties (Dependent)
        wave;           % Wavelength samples
        
        rows;           % number of rows in the cone mosaic
        cols;           % number of cols in the cone mosaic
        mosaicSize;     % [rows, cols]
        patternSupport; % [X(:) Y(:)] axes
        width;          % width of cone mosaic in meters
        height;         % height of cone mosaic in meters
        fov;            % horizontal/vertical field of view assuming inf 
                        % scene distance and 17mm optics focal length
        
        coneLocs;       % cone locations in meters
        qe;             % absorptance with macular pigment (not lens)
        
        current;        % The spatial array of photocurrent over time
        spatialDensity; % spatial density (ratio) of the K-LMS cones
    end
    
    properties (Access=private)
        % spatial density (ratio) of the K-LMS cones
        %
        % There are two properties about this LMS ratio: sptailDensity_ and
        % spatialDensity. spatialDensity_ is a private variable of the
        % class while spatialDensity is a dependent variable.
        %
        % When setting the spatialDensity_ inside the class, we just change
        % this private property directly. When setting spatialDensity from
        % outside, the set.spatialDensity is executed and the cone mosaic
        % is regenerated.
        %
        % It's possible to have set function for spatialDensity directly.
        % However, the current Matlab does not allow update other
        % properties (e.g. mosaic pattern) in the set function for
        % non-dependent properties.
        %
        % HJ
        spatialDensity_;
    end
    
    methods    
        % Constructor
        function obj = coneMosaic(varargin)
            % Initialize the cone mosaic class
            %   cMosaic =  coneMosaic('cone',cone,'os','os);
            
            % TODO:  We should take eccentricity and angle as an input
            % parameter and create cones of the appropriate size for that
            % retinal location
            
            % parse input
            p = inputParser;
            
            p.addParameter('name', 'cone mosaic', @ischar);
            
            % Included objects
            p.addParameter('pigment', photoPigment(), ...
                @(x) isa(x, 'photoPigment'));
            p.addParameter('macular', Macular(), @(x)isa(x, 'Macular'));
            p.addParameter('os', osLinear(), @(x)(isa(x,'outerSegment')));

            % Mosaic parameters
            p.addParameter('center',[0 0], @(x)(numel(x) ==2));
            p.addParameter('wave', 400:10:700, @isnumeric);
            p.addParameter('pattern', [], @isnumeric);
            p.addParameter('spatialDensity', [0 0.6 0.3 0.1], @isnumeric);
            p.addParameter('size', [72 88], @isnumeric);
            p.addParameter('integrationTime', 0.05, @isscalar);
            % p.addParameter('sampleTime', 0.001, @isscalar);
            
            % Computational features
            p.addParameter('emPositions', [0 0], @isnumeric);
            p.addParameter('noiseFlag', 1, @isscalar);
            
            p.parse(varargin{:});
            
            % set properties
            obj.name    = p.Results.name;
            obj.pigment = p.Results.pigment;
            obj.macular = p.Results.macular;
            obj.os      = p.Results.os;
            
            obj.center = p.Results.center(:)';
            obj.wave   = p.Results.wave;
            obj.spatialDensity_ = p.Results.spatialDensity(:);
            % obj.sampleTime      = p.Results.sampleTime;
            obj.integrationTime = p.Results.integrationTime;
            
            obj.noiseFlag = p.Results.noiseFlag;
            obj.emPositions = p.Results.emPositions;
            
            % Set the cone spacing and aperture given its eccentricity and
            % angle.  We could specify eye, but are we really sure about
            % the left right thing in human?
            % 
            % Units of returns are meters
            ecc = sqrt(sum(obj.center.^2));
            ang = atan2d(obj.center(2),obj.center(1));
            [spacing, aperture] = coneSize(ecc,ang);

            obj.pigment.pdWidth  = aperture;
            obj.pigment.pdHeight = aperture;
            obj.pigment.height = spacing;
            obj.pigment.width  = spacing;
            
            
            % Not sure why we do this.  Related to NP.  To discuss in an
            % issue. (BW)
            % From NP to BW: 
            % In a rect mosaic (instantiated via coneMosaic), obj.pattern is the actual pattern of
            % the cones in the mosaic. Therefore the elements of obj.pattern are separated by a 
            % distance equal to the inter-cone separation.
            % In a hex mosaic (instantiated via coneMosaicHex), obj.pattern represents the pattern 
            % of a high resolution sampling array, in which most elements are null(K) cones. 
            % The non-null elements are LMS cones positioned at the nodes of a hex grid. Therefore, 
            % in a coneMosaicHex, the patternSampleSize is set (by the coneMosaicHex class) to a 
            % value that corresponds not to actual LMS cone spacing, but to the spacing of the 
            % high resolution sampling array.
            obj.patternSampleSize = [obj.pigment.width, obj.pigment.height];
            
            % generate human cone mosaic pattern if not specified
            if isempty(p.Results.pattern)
                [~, obj.pattern] = humanConeMosaic(p.Results.size, ...
                    obj.spatialDensity_, obj.patternSampleSize(1));
            else
                obj.pattern = p.Results.pattern;
            end
            
            
            % Initialize the mosaic properties
            % obj.os.timeStep = obj.sampleTime;
            obj.os.patchSize = obj.width;
            
            % initialize listener
            % these listeners make sure the wavelength samples
            % in obj.pigment and obj.macular are the same
            addlistener(obj.pigment, 'wave', 'PostSet', @obj.setWave);
            addlistener(obj.macular, 'wave', 'PostSet', @obj.setWave);
            
            % Trying to remove sampleTime
            %  addlistener(obj, 'sampleTime', 'PostSet', @obj.setSampleTime);
            % addlistener(obj.os, 'timeStep', 'PostSet', @obj.setSampleTime);

        end
        
        
        
        clearData(obj, varargin);
        a = coneAbsorptions(obj,varargin);
        window(obj, varargin);
        obj = setSizeToFOV(obj, fov, varargin);
        
        %% get methods for dependent variables
        % http://www.mathworks.com/help/matlab/matlab_oop/specifying-methods-and-functions.html#bu4wzba
        % All functions that use dots in their names must be defined in the
        % classdef file, including:
        %   -  Converter methods that must use the package name as part of
        %      the class name because the class is contained in packages 
        %   -  Property set and get access methods
        function val = get.wave(obj)
            val = obj.pigment.wave;
        end
        
        function val = get.rows(obj)  % number of rows
            val = size(obj.pattern, 1);
        end
        
        function val = get.cols(obj)  % number of cols
            val = size(obj.pattern, 2);
        end
        
        function val = get.mosaicSize(obj)
            val = [obj.rows obj.cols];
        end
        
        function val = get.patternSupport(obj)
            x = (1:obj.cols) * obj.patternSampleSize(1); x = x - mean(x);
            y = (1:obj.rows) * obj.patternSampleSize(2); y = y - mean(y);
            [xx, yy] = meshgrid(x,y);
            val(:,:,1) = xx;
            val(:,:,2) = yy;
        end
        
        function val = get.width(obj)  % width of cone mosaic in meters
            val = obj.cols * obj.patternSampleSize(1);
        end
        
        function val = get.height(obj)  % height of cone mosaic in meters
            val = obj.rows * obj.patternSampleSize(2);
        end
        
        function val = get.fov(obj)  % horizontal/vertical field of view
            val = 2 * atand([obj.width obj.height]/2/0.017);
        end
        
        function val = get.coneLocs(obj) % cone locations in meters
            %
            x = (1:obj.cols) * obj.pigment.width; x = x - mean(x) + obj.center(1);
            y = (1:obj.rows) * obj.pigment.height; y = y - mean(y) + obj.center(2);
            
            [X, Y] = meshgrid(x, y);
            val = [X(:) Y(:)];
        end
        
        function val = get.qe(obj)
            % cMosaic.qe
            % compute effective absorptance with macular pigments
            %
            % The cone mosaic quantum efficiency is the product of the cone
            % photopigment absorptance times the macular pigment
            % transmittance times the cone photopigment peak absorption.
            %
            % The eye quantum efficiency, called in plot, includes the lens
            % transmittance.
            val = bsxfun(@times, obj.pigment.absorptance, ...
               obj.macular.transmittance)*diag(obj.pigment.peakEfficiency);
        end
        
        function val = get.spatialDensity(obj)
            val = obj.spatialDensity_;
        end
        
        function val = get.absorptions(obj)           
            val = double(obj.absorptions);
        end
        
        function val = get.current(obj)
            val = double(obj.os.coneCurrentSignal);
        end
        
        %% set method for class properties
        function set.spatialDensity(obj, val)
            if all(obj.spatialDensity_(:) == val(:)), return; end
            obj.spatialDensity_ = val;
            [~, obj.pattern] = humanConeMosaic(obj.mosaicSize, ...
                    val, obj.patternSampleSize(1));
            obj.clearData();
        end
        
        function set.wave(obj, val)
            obj.pigment.wave = val(:);
            obj.macular.wave = val(:);
        end
       
        function set.absorptions(obj, val)
            obj.absorptions = single(val);
        end
        
        function set.current(obj, val)
            obj.os.osSet('cone current signal', single(val));
        end
        
        function set.fov(obj, val) % set field of view
            obj.setSizeToFOV(val);
        end
        
        function set.mosaicSize(obj, val)
            if any(val ~= obj.mosaicSize)
                [~, obj.pattern] = humanConeMosaic(val, ...
                    obj.spatialDensity_, obj.patternSampleSize(1));
                obj.os.patchSize = obj.width;
                obj.clearData();
            end
        end
        
        function set.rows(obj, val)
            obj.mosaicSize = [val obj.cols];
        end
        
        function set.cols(obj, val)
            obj.mosaicSize = [obj.rows val];
        end
    end
    
    methods (Access=public)
        % Declare the compute method
        [absorptions, current] = compute(obj, oi, varargin);
        
        % Declare the computeCurrent method
        computeCurrent(obj, varargin);
    end

    methods (Static)
        [noisyImage, theNoise] = photonNoise(absorptions,varargin);
    end

    % Methods may be called by the subclasses, but are otherwise private
    methods (Access = protected)
        cpObj = copyElement(obj);
    end
    
    methods (Access = public, Hidden)
        absorptions = computeSingleFrame(obj, oi, varargin);
        absorptions = applyEMPath(obj, LMS, varargin);
    end        
    
    methods (Access = private)
        setWave(obj, src, ~)
    end
        
        %         function setSampleTime(obj, src, ~)
        %             warning('Sample time changed...can be bad');
        %             switch src.DefiningClass.Name
        %                 case 'coneMosaic'
        %                     obj.os.timeStep = obj.sampleTime;
        %                 otherwise
        %                     obj.sampleTime = obj.os.timeStep;
        %             end
        %         end
        %     end
end