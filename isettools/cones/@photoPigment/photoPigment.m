classdef photoPigment < hiddenHandle
    %photoPigment  Class for single cone photopigment and related properties
    %
    % Syntax:
    %    pigment = photoPigment;
    %
    % Description:
    %    This class contains properties for the photopigment absorption
    %    properties of a single cone cell. 
    %
    %    For the full cone mosaic, see the coneMosaic class
    %
    %    Most of the terms represented here are descriptions of the photopigment
    %    itself.  In addition, there are a few terms the capture the effective
    %    optical size of the photopigment absorption.
    %
    %    Default parameters are determined by underlying routines that get 
    %    the required data types. Unmatched key/value pairs passed to
    %    photoPigment are passed on to the underlying routines and can be
    %    used to adjust the parameters obtained. 
    %      absorbance    10.^getLogConeAborbance
    %
    %    [DHB NOTE: Need to explain about width and heigh, pdWidth and pdHeight and
    %    how these are used.  Perhaps even simplify code not to have both.]
    %
    % Input:
    %    None.
    %
    % Output:
    %    pigment           The created photoPigment object.
    %   
    % Optional key/value pairs:
    %
    %    'wave'              Vector of wavelengths in nm (400:10:31).
    %    'opticalDensity'    Three vector of optical densities for L, M and S cone photopigment (default: [0.5 0.5 0.4]).
    %    'absorbance'        L, M and S cone absorbance spectra. (Default empty, in which case these
    %                          are obtained through routine coneAbsorbanceReadData.)
    %    'peakEfficiency'    Peak quantal efficiency for isomerizations for L, M and S cones (default [2 2 2]/3).    
    %    'width'             Cone width (including gap between cones) in meters (default 2e-6).
    %    'height'            Cone height (including gap between cones) in meters  (default 2e-6).
    %    'pdWidth'           Collecting area width in meters (default 2e-6).
    %    'pdHeight'          Collecting area height in meters (default 2e-6).   
    %
    % See also: t_conePhotoPigment, coneMosaic, Macular, lens
    
    % HJ, ISETBIO Team, 2016
    
    properties  % public properties
        opticalDensity;  % photopigment optical densities for L,M,S
        peakEfficiency;  % peak absorptance efficiency
        
        width;           % cone width (including gap) in meters
        height;          % cone height (including gap) in meters
        
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
            p.KeepUnmatched = true;
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
            obj.wave_ = (390:830)';
            obj.opticalDensity = p.Results.opticalDensity(:);
            obj.peakEfficiency = p.Results.peakEfficiency(:);
            
            obj.width = p.Results.width;
            obj.height = p.Results.height;
            obj.pdWidth = p.Results.pdWidth;
            obj.pdHeight = p.Results.pdHeight;
            
            % If absorbance is not specified, we obtain it using the defaults
            % of coneAbsorbanceReadData.  
            if isempty(p.Results.absorbance)
                obj.absorbance_ = coneAbsorbanceReadData(p.Unmatched,'wave',obj.wave_);
            else
                obj.absorbance = p.Results.absorbance;
            end
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
        % When we have them, they go here
    end
end