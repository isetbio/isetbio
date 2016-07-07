classdef Lens
    % Class for human lens pigment properties
    %
    %     lens = Lens();
    %
    % The returned class object includes a variety of derived terms. This
    % should help to keep the relationship between entities straight.
    %
    % Useful formulae
    %
    %   Absorbance spectra are normalized to a peak value of 1.
    %   Absorptance spectra are the proportion of quanta actually absorbed.
    %   Equation: absorptanceSpectra = 1 - 10.^(-OD * absorbanceSpectra)
    %
    % The original lens densities values were taken from PTB and/or the
    % Stockman site. Go to http://cvision.ucsd.edu, then click on
    % Prereceptoral filters.
    %
    % Examples:
    %   lens = Lens();
    %
    % HJ/BW ISETBIO Team 2013.
    
    properties       % public properties
        name;        % Name of this particular lens object
        density;     % macular pigment density
        wave;        % wavelength samples in nm
    end
    
    properties (Access=private)
        wave_;         % internal wavelength samples
        unitDensity_;  % unit density absorbance sampled with wave_
    end
    
    properties (Dependent)
        unitDensity;     % spectral absorbance with unit pigment density
        spectralDensity; % unitDensity scaled by obj.density
        transmittance;   % proportion of quanta transmitted
        absorptance;     % proportion of quanta absorbed
    end
    
    methods  % public methods
        % constructor
        function obj = Lens(varargin)
            % thisLens = 
            %   Lens('wave',400:10:700, 'density',1,'name','my lens');
            %
            % parse input
            p = inputParser;
            p.addParameter('wave', 400:10:700, @isnumeric);
            p.addParameter('density', 1, @isscalar);
            p.addParameter('unitDensity', [], @isnumeric);
            p.addParameter('name','human lens',@isstr);
            
            p.parse(varargin{:});
            
            obj.name = p.Results.name;
            
            % set properties
            obj.wave = p.Results.wave(:);
            obj.wave_ = (400:10:700)';
            obj.density = p.Results.density;
            
            if isempty(p.Results.unitDensity)
                obj.unitDensity_ = ieReadSpectra('lensDensity.mat', ...
                    obj.wave_);
            else
                obj.unitDensity = p.Results.unitDensity;
            end
        end
        
        % get methods for dependent variables
        function val = get.unitDensity(obj)
            % interpolate for wavelength samples
            val = interp1(obj.wave_, obj.unitDensity_,obj.wave, 'pchip');
            val = max(val, 0);
        end
        
        function val = get.spectralDensity(obj)
            % compute scaled absorbance
            val = obj.unitDensity * obj.density;
        end
        
        function val = get.transmittance(obj)
            % compute proportion of quanta transmitted
            val = 10.^(-obj.spectralDensity);
        end
        
        function val = get.absorptance(obj)
            % comptue proportion of quanta absorbed
            val = 1 - obj.transmittance;
        end
        
        % set methods for dependent variables
        function obj = set.unitDensity(obj, val)
            % interpolate for wavelength samples
            obj.unitDensity_ = interp1(obj.wave, val, obj.wave_, 'pchip');
            obj.unitDensity_ = max(obj.unitDensity, 0);
        end
    end
end