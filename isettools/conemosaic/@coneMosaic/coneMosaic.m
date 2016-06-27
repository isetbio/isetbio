classdef coneMosaic < hiddenHandle
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
        
        pigment;          % Cone class object, contain single cone property
        macular;          % Macular class object
        os;               % Outersegment properties
        
        pattern;          % Pattern of K-LMS cones in the mosaick
        spatialDensity;   % spatial density (ratio) of the K-LMS cones
        integrationTime;  % Cone temporal integration time in secs
        sampleTime;       % Time step for em and os computation, shall we
                          % set this the same as integrationTime?
        
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
        fov;        % horizontal/vertical field of view assuming inf scene 
                    % distance and 17mm optics focal length
        
        coneLocs;   % cone locations in meters
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
            
            p.addParameter('pigment', photoPigment(), ...
                @(x) isa(x, 'photoPigment'));
            p.addParameter('macular', Macular(), @(x)isa(x, 'Macular'));
            p.addParameter('wave', 400:10:700, @isnumeric);
            p.addParameter('integrationTime', 0.05, @isscalar);
            p.addParameter('emPositions', [0 0], @isnumeric);
            p.addParameter('spatialDensity', [0 0.6 0.3 0.1], @isnumeric);
            p.addParameter('noiseFlag', 1, @isscalar);
            p.addParameter('pattern', [], @isnumeric);
            p.addParameter('size', [72 88], @isnumeric);
            p.addParameter('sampleTime', 0.001, @isscalar);
            p.addParameter('os', osCreate('linear'), ...
                @(x)(isa(x,'outerSegment')));
            
            p.parse(varargin{:});
            
            % set properties
            obj.pigment = p.Results.pigment;
            obj.macular = p.Results.macular;
            obj.os = p.Results.os;
            obj.wave = p.Results.wave;
            obj.integrationTime = p.Results.integrationTime;
            obj.spatialDensity = p.Results.spatialDensity(:);
            obj.noiseFlag = p.Results.noiseFlag;
            obj.sampleTime = p.Results.sampleTime;
            obj.emPositions = p.Results.emPositions;
            
            if isempty(p.Results.pattern)
                [~, obj.pattern] = humanConeMosaic(p.Results.size, ...
                    obj.spatialDensity, obj.pigment.width);
            else
                obj.pattern = p.Results.pattern;
            end
            
            % Initialize the mosaic properties
            % obj.matchSensor(varargin{:});
            
            % initialize listener
            % these listeners make sure the wavelength samples
            % in obj.pigment and obj.macular are the same
            addlistener(obj.pigment, 'wave', 'PostSet', @obj.setWave);
            addlistener(obj.macular, 'wave', 'PostSet', @obj.setWave);
        end
        
        % set size to fov
        function obj = setSizeToFOV(obj, fov, varargin)
            % set cone mosaic size according to the field of view
            %
            % Inputs:
            %   fov    - 2-element vector for desired horizontal/vertical
            %            field of view in degrees
            % 
            % Outputs:
            %   obj    - class object with size and mosaic updated
            %
            
            % parse input
            p = inputParser;
            p.addRequired('fov', @isnumeric);
            p.addParameter('sceneDist', inf, @isnumeric);
            p.addParameter('focalLength', 0.017, @isnumeric);
            
            p.parse(fov, varargin{:});
            sceneDist = p.Results.sceneDist;
            focalLength = p.Results.focalLength;
            if isscalar(fov), fov = [fov fov]; end
            fov = fov([2 1]); % flip the order to be vertical/horizontal
            
            % compute current field of view
            imageDist = 1/(1/focalLength - 1/sceneDist);
            curFov = 2*atand([obj.height obj.width]/imageDist/2);
            
            % set new size to object
            obj.mosaicSize = ceil(obj.mosaicSize .* fov./curFov);
        end
        
        
        % get methods for dependent variables
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
        
        function val = get.width(obj)  % width of cone mosaic in meters
            val = obj.cols * obj.pigment.width;
        end
        
        function val = get.height(obj)  % height of cone mosaic in meters
            val = obj.rows * obj.pigment.height;
        end
        
        function val = get.fov(obj)  % horizontal/vertical field of view
            val = 2 * atand([obj.width obj.height]/2/0.017);
        end
        
        function val = get.coneLocs(obj) % cone locations in meters
            x = (1:obj.cols) * obj.pigment.width; x = x - mean(x);
            y = (1:obj.rows) * obj.pigment.height; y = y - mean(y);
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
        
        % set method for dependent variables
        function set.wave(obj, val)
            obj.pigment.wave = val(:);
            obj.macular.wave = val(:);
        end
        
        function set.fov(obj, val) % set field of view
            obj.setSizeToFOV(val);
        end
        
        function set.mosaicSize(obj, val)
            if any(val ~= obj.mosaicSize)
                [~, obj.pattern] = humanConeMosaic(val, ...
                    obj.spatialDensity, obj.pigment.width);
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
        
        function [absorptions, current] = compute(obj, oi, varargin)
            % coneMosaic.plot()
            %
            % Compute the pattern of cone absorptions and typically the
            % photocurrent.
            p = inputParser;
            p.addRequired('oi',@isstruct);
            p.addParameter('currentFlag', false, @islogical);
            
            % newNoise false means frozen noise, using the seed rng(1)
            p.addParameter('newNoise',true,@islogical);
            
            p.parse(oi,varargin{:});
            oi = p.Results.oi;
            currentFlag = p.Results.currentFlag;
            newNoise = p.Results.newNoise;
            
            % Deal with eye movements
            % prepare parameters for eye movement
            if isempty(obj.emPositions), obj.emPositions = [0 0]; end
            mask = zeros(obj.rows, obj.cols, 3); % locations for cones
            for ii = 2 : 4 % L, M, S
                mask(:,:,ii-1) = double(obj.pattern == ii);
            end
            
            % extend sensor size
            xpos = obj.emPositions(:, 1); ypos = obj.emPositions(:, 2);
            padRows = max(abs(ypos)); padCols = max(abs(xpos));
            cpObj = obj.copy();
            cpObj.pattern = zeros(obj.rows+2*padRows, obj.cols+2*padCols);
            
            % compute full LMS noise free absorptions
            LMS = cpObj.computeSingleFrame(oi, 'fullLMS', true);
            absorptions = zeros(obj.rows,obj.cols,size(obj.emPositions,1));
            
            for ii = 1 : size(obj.emPositions, 1)
                % select out the subframe given the eye position
                cropLMS = LMS(1+padRows+ypos(ii):end-padRows+ypos(ii), ...
                    1+padCols-xpos(ii):end-padCols-xpos(ii), :);
                
                % sample by conetype
                absorptions(:, :, ii) = sum(cropLMS .* mask, 3);
            end
            
            % Add photon noise to the whole volume
            if obj.noiseFlag
                absorptions = obj.photonNoise(absorptions, ...
                    'newNoise', newNoise);
            end
            obj.absorptions = absorptions;
            
            % If we want the photo current, use the os model
            if currentFlag
                current = obj.os.compute;
                obj.current = current;
            end
        end
    end
    
    methods (Static)
        function [noisyImage, theNoise] = photonNoise(absorptions,varargin)
            % Photon noise at the absorptions is Poisson.
            % The Poisson variance is equal to the mean
            % For moderate mean levels, Poisson is very close to Normal
            %
            % We multiply each point in the absorption data by the square
            % root of its mean value to create the noise standard
            % deviation. For most cases the Normal approximation is
            % adequate. We trap (below) the cases when the value is small
            % (less than 25) and replace it with the real Poisson random
            % value, which is slower to compute.
            
            p = inputParser;
            p.addRequired('absorptions',@isnumeric);
            
            % Frozen or new random noise
            p.addParameter('newNoise',true,@islogical);
            
            p.parse(absorptions,varargin{:});
            absorptions  = p.Results.absorptions;
            newNoise = p.Results.newNoise;
            
            % This is std * N(0,1)
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
        
        function pos = emGenSequence(nFrames, varargin)
            % Generate eye movement path
            %
            % Inputs
            %   nFrames - number of frames to generate
            %
            % Optional inputs (key value pairs in varargin):
            %   'em'    - eye movement structure, see emCreate for details
            %   'rSeed' - random seed to be used
            %
            % Ouputs
            %   pos     - nFramesx2 matrix of eye positions in units of
            %             cone positions
            %
            % See also:
            %   emCreate
            
            % parse input
            p = inputParser;
            p.addRequired('nFrames', @isscalar);
            p.addParameter('em', emCreate, @isstruct);
            p.addParameter('rSeed', [], @isscalar);
            
            % set parameters
            p.parse(nFrames, varargin{:});
            em = p.Results.em;
            if ~isempty(p.Results.rSeed), rng(p.Results.rSeed); end
            
            emFlag = emGet(em, 'em flag');
            pos = zeros(nFrames, 2);
            sampTime  = emGet(em, 'sample time');
            
            
            % generate eye movement for tremor
            if emFlag(1)
                % Load parameters
                amplitude  = emGet(em, 'tremor amplitude', 'cones/sample');
                interval   = emGet(em, 'tremor interval');
                intervalSD = emGet(em, 'tremor interval SD');
                
                % Compute time of tremor occurs
                t = interval + randn(nFrames, 1) * intervalSD;
                t(t < 0.001) = 0.001; % get rid of negative values
                tPos = cumsum(t);
                tPos = round(tPos / sampTime);
                indx = 1:find(tPos <= nFrames, 1, 'last');
                tPos = tPos(indx);
                
                % Generate random step at the selected times
                direction = rand(length(tPos),1);
                
                % Unit length direction
                pos(tPos, :) = amplitude*[direction sqrt(1-direction.^2)];
                
                pos(tPos, :) = bsxfun(@times,pos(tPos,:),t(indx)/sampTime);
                pos = pos .* (2*(randn(size(pos))>0)-1); % shuffle the sign
                pos = cumsum(pos, 1);
            end
            
            % generate eye movement for drift
            if emFlag(2)
                % Load Parameters
                speed     = emGet(em, 'drift speed', 'cones/sample');
                speedSD   = emGet(em, 'drift speed SD', 'cones/sample');
                
                % Generate random move at each sample time
                theta = 360 * randn + 0.1 * (1 : nFrames)';
                direction = [cosd(theta) sind(theta)];
                s = speed + speedSD * randn(nFrames, 1);
                pos = filter(1,[1 -1],bsxfun(@times, direction, s)) + pos;
            end
            
            % generate eye movement for micro-saccade
            if emFlag(3)
                % Load parameters
                interval = emGet(em, 'msaccade interval');
                intervalSD = emGet(em, 'msaccade interval SD');
                dirSD = emGet(em, 'msaccade dir SD', 'deg');
                speed = emGet(em, 'msaccade speed', 'cones/sample');
                speedSD = emGet(em, 'msaccade speed SD', 'cones/sample');
                
                % compute time of occurance
                t = interval + randn(nFrames, 1) * intervalSD;
                t(t < 0.3) = 0.3 + 0.1*rand; % get rid of negative values
                t = cumsum(t);
                tPos = round(t / sampTime);
                tPos = tPos(1:find(tPos <= nFrames, 1, 'last'));
                
                % Compute positions
                for ii = 1 : length(tPos)
                    curPos = pos(tPos(ii), :);
                    duration = round(sqrt(curPos(1)^2+curPos(2)^2)/speed);
                    direction = atand(curPos(2)/curPos(1)) + dirSD * randn;
                    direction = [cosd(direction) sind(direction)];
                    direction = abs(direction) .* (2*(curPos < 0) - 1);
                    
                    offset = zeros(nFrames, 2);
                    indx = tPos(ii):min(tPos(ii) + duration - 1, nFrames);
                    curSpeed = speed + speedSD * randn;
                    if curSpeed < 0, curSpeed = speed; end
                    offset(indx, 1) = curSpeed*direction(1);
                    offset(indx, 2) = curSpeed*direction(2);
                    
                    pos = pos + cumsum(offset);
                end
            end
            pos = round(pos);           
        end
    end
    
    % Methods may be called by the subclasses, but are otherwise private
    methods (Access = protected)
        function cpObj = copyElement(obj)
            % make a shallow copy of the properties
            cpObj = copyElement@matlab.mixin.Copyable(obj);
            
            % make deep copy of the cone and macular class
            cpObj.pigment = cpObj.pigment.copy();
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
            p.addParameter('fullLMS', false, @islogical); % full LMS images
            
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
                oiGet(oi, 'height spatial resolution'), ...
                oiGet(oi, 'width spatial resolution'));
            [coneR, coneC] = sample2space(0.5:obj.rows - 0.5, ...
                0.5:obj.cols - 0.5, obj.pigment.height, obj.pigment.width);
            
            if fullLMS
                absDensity = zeros(obj.rows, obj.cols, 3);
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
            absorptions=absDensity*obj.pigment.pdArea*obj.integrationTime;
        end
        
        % callback function for listeners
        % When you set the cone, it sets the macular, and vice versa.
        function setWave(obj, src, ~)
            switch src.DefiningClass.Name
                case 'photoPigment'
                    obj.macular.wave = obj.pigment.wave;
                case 'Macular'
                    obj.pigment.wave = obj.macular.wave;
            end
        end
    end
    
end