classdef coneMosaic < handle & matlab.mixin.Copyable
    % Create a cone mosaic class
    %
    %   cMosaic =  coneMosaic('cone', cone, 'os', os);
    %
    % The cone mosaic defines an array of cones.  The individual cones have
    % absorption properties defined by cMosaic.cone.  The computation from
    % absorptions to photocurrent is defined by cMosaic.os
    %
    % HJ/JRG/BW ISETBIO Team, 2016
    
    properties  % public properties
        
        cone;             % Cone class object, contain single cone property
        macular;          % Macular class object
        os;               % Outersegment properties
        
        pattern;          % Pattern of K-LMS cones in the mosaick
        integrationTime;  % Cone temporal integration time in secs
        emPositions;      % Eye movement positions in number of cones.
                          % The length of this property controls number of
                          % frames to be computed
        noiseFlag;        % To control which noise is included
    end
    
    properties (GetAccess=public, SetAccess=private)
        absorptions;      % The spatial array of cone absorptions over time
        current;          % The spatial array of photocurrent over time
    end
    
    properties (Dependent)
        wave;       % Wavelength samples
        
        rows;       % number of rows in the cone mosaic
        cols;       % number of cols in the cone mosaic
        mosaicSize; % [rows, cols]
        
        width;      % width of cone mosaic in meters
        height;     % height of cone mosaic in meters
        
        qe;         % effective absorptance with macular pigment (not lens)
    end
    
    % Public methods
    methods
        
        % Constructor
        function obj = coneMosaic(varargin)
            % Initialize the cone mosaic class
            %
            %   cMosaic =  coneMosaic('cone',cone,'os','os);
            %
            
            % parse input
            p = inputParser;
            
            p.addParameter('cone', Cones(),@(x) isa(x, 'Cones'));
            p.addParameter('macular', Macular(), @(x)isa(x, 'Macular'));
            p.addParameter('wave', 400:10:700, @isnumeric);
            p.addParameter('integrationTime', 0.05, @isscalar);
            p.addParameter('emPositions', [0 0], @isnumeric);
            p.addParameter('noiseFlag', 1, @isscalar);
            p.addParameter('pattern', [], @isnumeric);
            p.addParameter('size', [72 88], @isnumeric);
            p.addParameter('os', osCreate('linear'), ...
                @(x)(isa(x,'outerSegment')));
            
            p.parse(varargin{:});
            
            % set properties
            obj.cone = p.Results.cone;
            obj.macular = p.Results.macular;
            obj.os = p.Results.os;
            obj.wave = p.Results.wave;
            obj.integrationTime = p.Results.integrationTime;
            obj.noiseFlag = p.Results.noiseFlag;
            obj.emPositions = p.Results.emPositions;
            
            if isempty(p.Results.pattern)
                [~, obj.pattern] = humanConeMosaic(p.Results.size, ...
                    obj.cone.spatialDensity, obj.cone.width);
            else
                obj.pattern = p.Results.pattern;
            end
            
            % Initialize the mosaic properties
            % obj.matchSensor(varargin{:});
            
            % initialize listener
            % these listeners make sure the wavelength samples
            % in obj.cone and obj.macular are the same
            addlistener(obj.cone, 'wave', 'PostSet', @obj.setWave);
            addlistener(obj.macular, 'wave', 'PostSet', @obj.setWave);
        end
        
        
        % get methods for dependent variables
        function val = get.wave(obj)
            val = obj.cone.wave;
        end
        
        function val = get.rows(obj)  % number of rows
            val = size(obj.pattern, 1);
        end
        
        function val = get.cols(obj)  % number of cols
            val = size(obj.pattern, 2);
        end
        
        function val = get.width(obj)  % width of cone mosaic in meters
            val = obj.cols * obj.cone.width;
        end
        
        function val = get.height(obj)  % height of cone mosaic in meters
            val = obj.rows * obj.cone.height;
        end
        
        function val = get.qe(obj)
            % compute effective absorptance with macular pigments
            val = bsxfun(@times, obj.cone.absorptance, ...
                obj.macular.transmittance) * diag(obj.cone.peakEfficiency);
        end
        
        % set method for dependent variables
        function set.wave(obj, val)
            obj.cone.wave = val(:);
            obj.macular.wave = val(:);
        end
        
        function set.mosaicSize(obj, val)
            if any(val ~= obj.mosaicSize)
                [~, obj.pattern] = humanConeMosaic(val, ...
                    obj.cone.spatialDensity, obj.cone.width); 
            end
        end
        
        function set.rows(obj, val)
            obj.mosaicSize = [val obj.cols];
        end
        
        function set.cols(obj, val)
            obj.mosaicSize = [obj.rows val];
        end
        
        % set function, see osLinearSet for details
        % function obj = set(obj, varargin)
        % end
        
        % get function, see osLinearGet for details
        % function val = get(obj, varargin)
        %    val = obj;  % Place holder for what will happen!
        % end
        
    end
    
    % Methods that must only be implemented (Abstract in parent class).
    methods (Access=public)
        
        function obj = compute(obj, oi, varargin)
            % coneMosaic.plot()
            %
            % Compute the pattern of cone absorptions and typically the
            % photocurrent.  If you don't want the current, say so.
            p = inputParser;
            p.addRequired('oi',@isstruct);
            p.addParameter('current', false, @islogical); % Flag to compute current
            
            % newNoise false means frozen noise, using the seed rng(1)
            p.addParameter('newNoise',true,@islogical);   % Controls photon noise
            
            p.parse(oi,varargin{:});
            oi = p.Results.oi;
            currentFlag = p.Results.current;
            newNoise = p.Results.newNoise;
            
            % Deal with eye movements
            if isempty(obj.emPositions), obj.emPositions = [0 0]; end
            
            % First compute a larger frame of the LMS cones that extends as
            % far as we will need given the eye movement
            obj.extendConeArray;
            
            % Noise free frame
            obj.absorptions = obj.computeSingleFrame(oi);
            
            for ii=1:size(obj.emPositions,1)
                % Select out the subframe given the eye position
                obj.absorptions = obj.eyeposition;
            end
            
            % Add photon noise to the whole volume
            if obj.noiseFlag
                obj.absorptions = obj.photonNoise(obj.absorptions,'newNoise',newNoise);
            end
            
            % If we want the photo current, use the os model
            if currentFlag, obj.os.compute; end
        
        end
        
        function plot(obj, type, varargin)
            % coneMosaic.plot()
            % Do some plotting based on the input arguments.
        end
        
        function [noisyImage, theNoise] = photonNoise(obj,absorptions,varargin)
            % Photon noise at the absorptions is Poisson.
            % The Poisson variance is equal to the mean
            % For moderate mean levels, Poisson is very close to Normal
            % Randn is unit normal (N(0,1)).
            % S*Randn is N(0,S).
            %
            % We multiply each point in the absorption data by the square
            % root of its mean value to create the noise standard
            % deviation. For most cases the Normal approximation is
            % adequate. We trap (below) the cases when the value is small
            % (less than 25) and replace it with the real Poisson random
            % value, which is slower to compute.
            
            p = inputParser;
            p.addRequired('absorptions',@isnumeric);  
            p.addParameter('newNoise',true,@islogical);  % Frozen or new random noise
            
            p.parse(absorptions,varargin{:});
            absorptions  = p.Results.absorptions;
            newNoise = p.Results.newNoise;
            
            % This is S*N(0,1) = N(0,S)
            theNoise = sqrt(absorptions) .* randn(size(absorptions));
            
            % We add the mean electron and noise electrons together.
            noisyImage = round(absorptions + theNoise);
            
            % Now, we find the small mean values and create a true Poisson
            % sample. The Poisson algorithm is slow for big numbers, but it
            % is fast enough for small numbers. We can't rely on the Stats
            % toolbox being present, so we use this Poisson sampler from
            % Knuth. 
            
            % Apply the Poisson when the mean is less than this
            poissonCriterion = 25;
            idx = find(absorptions(:) < poissonCriterion);
            v = absorptions(absorptions < poissonCriterion);
            if ~isempty(v)
                vn = iePoisson(v, 1, newNoise);  % Poisson samples
                % for ii=1:length(r)
                %    theNoise(r(ii),c(ii))   = vn(ii);
                % For low mean values, we *replace* the mean value with the Poisson
                % noise; we do not *add* the Poisson noise to the mean
                %    noisyImage(r(ii),c(ii)) = vn(ii);
                %end
                theNoise(idx) = vn - absorptions(idx);
                noisyImage(idx) = vn;
            end
            
            % Find the highly unusual case in which the sum of the mean and
            % noise are less than zero. This only happens if the mean is 25
            % or greater and the noise is less than -25.  Very unusual, but
            % it can happen given how many times we run this random number
            % generate.  So this little clipping is an imperfection but it
            % isn't the worst thing we do.
            idx = (noisyImage < 0);
            noisyImage(idx)  = 0;
            theNoise(idx)    = noisyImage(idx) - absorptions(idx);
        end
        
    end
    
    % Methods may be called by the subclasses, but are otherwise private
    methods (Access = protected)
        function cpObj = copyElement(obj)
            % make a shallow copy of the properties
            cpObj = copyElement@matlab.mixin.Copyable(obj);
            
            % make deep copy of the cone and macular class
            cpObj.cone = cpObj.cone.copy();
            cpObj.macular = cpObj.macular.copy();
        end
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
        % compute function for single frame
        function absorptions = computeSingleFrame(obj, oi, varargin)
            % This function computes mean expected photon absorptions for
            % one frame, without including eye movements or noise.
            
            % parse inputs
            p = inputParser();
            p.addRequired('oi', @isstruct);
            p.addParameter('fullLMS', false, @islogical);  % No sampling, full LMS images
            
            p.parse(oi, varargin{:});
            
            fullLMS = p.Results.fullLMS;  % Logical
            
            % make a copy of current obj and set cone wavelength samples to
            % be same as oi
            obj = obj.copy();
            obj.wave = oiGet(oi, 'wave');
            
            % get scaled spectral qe, which includes cone pigment and
            % macular pigment properties. (Lens is in oi).
            sQE = obj.qe * oiGet(oi, 'bin width');
            
            % compute cone absorption density at oi sampled locations
            [photons, r, c] = RGB2XWFormat(oiGet(oi, 'photons'));
            absDensityLMS = XW2RGBFormat(photons * sQE, r, c);
            
            % regrid the density from oi sample locations to cone locations
            [oiR, oiC] = sample2space(0:r-1, 0:c-1, ...
                oiGet(oi, 'height spatial resolution'), oiGet(oi, 'width spatial resolution'));
            [coneR, coneC] = sample2space(0.5:obj.rows - 0.5, ...
                0.5:obj.cols - 0.5, obj.cone.height, obj.cone.width);
            
            if fullLMS
                absDensity = zeros(obj.rows, obj.height, 3);
            else
                absDensity = 0;
            end
            warning('off','MATLAB:interp1:NaNinY');
            for ii = 2 : 4  % loop through L, M and S, 1 = Blank/Black
                curDensity = interp1(oiR, absDensityLMS(:,:,ii-1), ...
                    coneR, 'linear', 0)';
                curDensity = interp1(oiC, curDensity, coneC, 'linear', 0)';
                if fullLMS
                    absDensity(:, :, ii-1) = curDensity;
                else
                    % Pick out the relevant cones by position
                    absDensity = absDensity+(obj.pattern==ii).*curDensity;
                end
            end
            warning('on','MATLAB:interp1:NaNinY');
            
            % Sometimes we don't have the cone type so we have a bad
            % number. Set the missing values to 0
            absDensity(isnan(absDensity)) = 0;
            
            % compute expected cone absorptions
            absorptions = absDensity*obj.cone.pdArea*obj.integrationTime;
        end
        
        % callback function for listeners
        % When you set the cone, it sets the macular, and vice versa.
        function setWave(obj, src, ~)
            switch src.DefiningClass.Name
                case 'Cones'
                    obj.macular.wave = obj.cone.wave;
                case 'Macular'
                    obj.cone.wave = obj.macular.wave;
            end
        end
    end
    
end