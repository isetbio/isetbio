classdef Lens < handle
    % LENS  Class for human lens pigment properties
    %     lens = LENS() returns a lens object with a variety of derived
    %     properties.
    %
    %     Useful factoids:
    %     Absorbance spectra are normalized to a peak value of 1.
    %     Absorptance spectra are the proportion of quanta actually absorbed.
    %     Equation: absorptanceSpectra = 1 - 10.^(-OD * absorbanceSpectra)
    %
    %     The original lens densities values were taken from PTB and/or the
    %     Stockman site. Go to http://cvision.ucsd.edu, then click on
    %     Prereceptoral filters.
    %
    %     Examples:
    %
    %         lens = Lens();
    
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
            obj.wave_ = (390:830)';
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
            % Can interpolate
            val = 10.^(-obj.spectralDensity);
        end
        
        function val = get.absorptance(obj)
            % compute proportion of quanta absorbed
            val = 1 - obj.transmittance;
        end
        
        function val = get.density(obj)
            val = obj.density;
        end
        
        function set(obj,param,val)
            p = inputParser; 
            p.KeepUnmatched = true;
            p.addRequired('param', @isstr);
            p.addRequired('val');
            
            % set parameter value
            switch ieParamFormat(param)
                case 'name'
                    assert(ischar(val),'Name should be a string');
                    obj.name = val;
                case {'wave', 'wavelength'}
                    assert(isvector(val), 'wave should be vector');
                    obj.wave = val(:);
                case {'absorbance','unitdensity'}
                    assert(length(val) == length(lens.wave), ...
                        'Val should have same length as lens wavelength');
                    obj.unitDensity_ = interp1(obj.wave, val, obj.wave_, 'pchip');
                    obj.unitDensity_ = max(obj.unitDensity, 0);
                case 'density'
                    assert(isscalar(val), 'val should be scalar');
                    obj.density = val;
                otherwise
                    error('Unknown parameter %s\n',param);
            end
        end
        
        % set methods for dependent variables
        %         function obj = set.unitDensity(obj, val)
        %             % interpolate for wavelength samples
        %             obj.unitDensity_ = interp1(obj.wave, val, obj.wave_, 'pchip');
        %             obj.unitDensity_ = max(obj.unitDensity, 0);
        %         end
        
        %         function set.wave(obj, val)
        %             obj.wave = val;
        %         end
        
        %         function set.density(obj, val)
        %             obj.density = val;
        %         end
        
    end
end