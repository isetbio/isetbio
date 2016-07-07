classdef photoPigment < hiddenHandle
    % Class for single cone
    %
    %   cone = Cones();
    %
    % This class contains properties for single cone cell. For the full
    % cone mosaic, see coneMosaic class
    %
    % HJ, ISETBIO Team, 2016
    
    properties  % public properties
        opticalDensity;  % photopigment optical densities for L,M,S
        peakEfficiency;  % peak absorptance efficiency
        width;           % cone width (include gap) in meters
        height;          % cone height (include gap) in meters
        
        pdWidth;         % photodetector width in meters
        pdHeight;        % photodetector height in meters
    end
    
    properties (SetObservable, AbortSet)
        wave;            % wavelength samples
    end
    
    properties (Dependent)
        absorbance;         % spectral absorbance of the cones
        absorptance;        % cone absorptance without ocular media
        quantaFundamentals; % normalized cone absorptance
        energyFundamentals; % normalized cone absorption in energy units
        area;               % width * height
        pdArea;             % pdWidth * pdHeight
        gapWidth;           % width - pdWidth
        gapHeight;          % height - pdHeight
    end
    
    properties(Access = private)  % private properties
        wave_;        % internal wavelength samples
        absorbance_;  % absorbance data sampled at wave_
    end
    
    methods  % public methods
        % constructor
        function obj = photoPigment(varargin)
            % parse input
            p = inputParser;
            p.addParameter('wave', 400:10:700, @isnumeric);
            p.addParameter('opticalDensity', [0.5 0.5 0.4], @isnumeric);
            p.addParameter('absorbance', [], @isnumeric);
            p.addParameter('peakEfficiency', [2 2 2]/3, @isnumeric);
            p.addParameter('width', 2e-6, @isnumeric);
            p.addParameter('height', 2e-6, @isnumeric);
            p.addParameter('pdWidth', 2e-6, @isnumeric);
            p.addParameter('pdHeight', 2e-6, @isnumeric);
            
            p.parse(varargin{:});
            
            % set object properties
            obj.wave = p.Results.wave(:);
            obj.wave_ = (400:10:700)';
            obj.opticalDensity = p.Results.opticalDensity(:);
            obj.peakEfficiency = p.Results.peakEfficiency(:);
            
            obj.width = p.Results.width;
            obj.height = p.Results.height;
            obj.pdWidth = p.Results.pdWidth;
            obj.pdHeight = p.Results.pdHeight;
            
            if isempty(p.Results.absorbance)
                obj.absorbance_ = 10 .^ ...
                    ieReadSpectra('coneAbsorbance', obj.wave_);
            else
                obj.absorbance = p.Results.absorbance;
            end
        end
        
        function str = description(obj, varargin)
            % generate description string for this object
            str = sprintf('\tWidth/height:\t[%.2f, %.2f] um\n', ...
                obj.pdWidth*1e6, obj.pdHeight*1e6);
            str = [str sprintf('\tGap (h, v):\t[%.2f, %.2f] um\n', ...
                obj.gapWidth*1e6, obj.gapHeight*1e6)];
            str = [str sprintf('\tOptical density: ') ...
                sprintf('[%.2f, %.2f, %.2f]\n', obj.opticalDensity(1), ...
                obj.opticalDensity(2), obj.opticalDensity(3))];
            str = [str sprintf('\tPeak efficiency: ') ...
                sprintf('[%.2f, %.2f, %.2f]\n', obj.peakEfficiency(1), ...
                obj.peakEfficiency(2), obj.peakEfficiency(3))];
        end
        
        % get method for dependent variable
        function val = get.absorbance(obj) % inerpolate for absorbance
            val = interp1(obj.wave_, obj.absorbance_, obj.wave, ...
                'linear', 'extrap');
            val = ieClip(val, 0, 1);
        end
        
        function val = get.absorptance(obj) % compute absorptance
            val = 1 - 10.^(-obj.absorbance*diag(obj.opticalDensity));
        end
        
        function val = get.quantaFundamentals(obj)
            % compute quanta fundamentals
            val = bsxfun(@rdivide, obj.absorptance, max(obj.absorptance));
        end
        
        function val = get.energyFundamentals(obj)
            h = vcConstants('planck');
            c = vcConstants('speed of light');
            val = 1e-9*bsxfun(@times,obj.quantaFundamentals/h/c,obj.wave);
        end
        
        function val = get.area(obj)
            val = obj.width * obj.height;
        end
        
        function val = get.gapWidth(obj)
            val = obj.width - obj.pdWidth;
        end
        
        function val = get.gapHeight(obj)
            val = obj.height - obj.pdHeight;
        end
        
        function val = get.pdArea(obj)
            val = obj.pdWidth * obj.pdHeight;
        end
        
        % set method for dependent variable
        function set.absorbance(obj, val)
            obj.absorbance_ = interp1(obj.wave, val, obj.wave_, ...
                'linear', 'extrap');
            obj.absorbance_ = ieClip(obj.absorbance_, 0, 1);
        end
    end
    
    methods (Static)
        function density = eccDensity(eccMM, angle, whichEye, varargin)
            % Compute cone packing density as a function of retinal
            % position
            %
            %   density = coneDensity(ecc, angle, whichEye, varargin)
            %
            % Inputs:
            %   ecc      - eccentricity (retinal position amplitude) in mm
            %   angle    - retinal position angle in degree, default is 0
            %   whichEye - can be either 'left' (default) or 'right'
            %
            % Outputs:
            %   density  - cone packing density in cones/mm^2
            %
            % References:
            %   1) Curcio, C. A., Sloan, K. R., Kalina, R. E. and
            %      Hendrickson, A. E. (1990), Human photoreceptor
            %      topography. J. Comp. Neurol., 292: 497?523. doi:
            %      10.1002/cne.902920402
            %   2) Song, H., Chui, T. Y. P., Zhong, Z., Elsner, A. E., &
            %      Burns, S. A. (2011). Variation of Cone Photoreceptor
            %      Packing Density with Retinal Eccentricity and Age.
            %      Investigative Ophthalmology & Visual Science, 52(10),
            %      7376?7384. http://doi.org/10.1167/iovs.11-7199
            %
            % Example:
            %   density = Cones.eccDensity(0, 90, 'left');
            %
            % HJ, ISETBIO TEAM, 2015
            
            % Check inputs
            if notDefined('eccMM'), eccMM = 0; end
            if notDefined('angle'), angle = 0; end
            if notDefined('whichEye'), whichEye = 'left'; end
            
            % load data
            d = load('coneDensity.mat');
            
            % interpolate for retinal position amplitude on axis (nasal,
            % superior, temporal and inferior direction)
            onAxisD = zeros(5, 1);
            angleQ = [0 90 180 270 360];
            
            % compute packing density for superior and inferior
            onAxisD(2)=interp1(d.superior.eccMM,d.superior.density,eccMM);
            onAxisD(4)=interp1(d.inferior.eccMM,d.inferior.density,eccMM);
            
            % nasal and temporal
            switch lower(whichEye)
                case 'left'
                    onAxisD(1) = interp1(d.nasal.eccMM, ...
                        d.nasal.density, eccMM);
                    onAxisD(3) = interp1(d.temporal.eccMM, ...
                        d.temporal.density, eccMM);
                case 'right'
                    onAxisD(1) = interp1(d.temporal.eccMM, ...
                        d.temporal.density, eccMM);
                    onAxisD(3) = interp1(d.nasal.eccMM, ...
                        d.nasal.density, eccMM);
                otherwise
                    error('unknown input for whichEye');
            end
            onAxisD(5) = onAxisD(1);
            
            % Interpolate for angle
            density = interp1(angleQ, onAxisD, angle, 'linear');
        end
    end
end