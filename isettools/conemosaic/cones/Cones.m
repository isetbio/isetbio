classdef Cones < hiddenHandle
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
        spatialDensity;  % spatial density (ratio) of the K-LMS cones
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
        function obj = Cones(varargin)
            % parse input
            p = inputParser;
            p.addParameter('wave', 400:10:700, @isnumeric);
            p.addParameter('opticalDensity', [0.5 0.5 0.4], @isnumeric);
            p.addParameter('absorbance', [], @isnumeric);
            p.addParameter('peakEfficiency', [2 2 2]/3, @isnumeric);
            p.addParameter('spatialDensity', [0 0.6 0.3 0.1], @isnumeric);
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
            obj.spatialDensity = p.Results.spatialDensity(:);
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
        
        % get method for dependent variable
        function val = get.absorbance(obj) % inerpolate for absorbance
            val = interp1(obj.wave_,obj.absorbance_,obj.wave,'linear',0);
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
            obj.absorbance_ = interp1(obj.wave,val,obj.wave_,'linear',0);
        end
    end
    
end