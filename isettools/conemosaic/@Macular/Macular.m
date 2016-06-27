classdef Macular < hiddenHandle
    % Class for macular pigment properties
    %
    %     macular = Macular('wave', wave)
    %
    % The human retina contains a pigment that covers the central (macular)
    % region. This macular pigment passes certain wavelengths of light more
    % than others. The pigment varies in density from central vision, where
    % it is highest, to increasingly peripheral vision.
    %
    % This function returns several measures of the macular pigment
    % wavelength properties as a function of macular pigment density (high
    % in the fovea, lower in the near fovea).
    %
    % The returned class object includes a variety of derived terms.
    % This should help to keep the relationship between entities straight.
    %
    % Density is the estimated (average) peak density of the pigment across
    % a variety of observers.  They estimate the average (across observers)
    % peak density to be 0.28, with a range of 0.17 to 0.48.
    %
    % Useful formulae
    %
    %   Absorbance spectra are normalized to a peak value of 1.
    %   Absorptance spectra are the proportion of quanta actually absorbed.
    %   Equation: absorptanceSpectra = 1 - 10.^(-OD * absorbanceSpectra)
    %
    % The original macular densities values were taken from the Stockman
    % site. Go to http://cvision.ucsd.edu, then click on Prereceptoral
    % filters.  At this point in time, I think the Psychtoolbox and the new
    % Stockman site are authoritative.
    %
    %
    % Examples:
    %   macular = Macular();
    %
    % HJ, ISETBIO TEAM, 2016
    
    properties       % public properties
        density;     % macular pigment density
    end
    
    properties (SetObservable, AbortSet)
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
        function obj = Macular(varargin)
            % parse input
            p = inputParser;
            p.addParameter('wave', 400:10:700, @isnumeric);
            p.addParameter('density', 0.35, @isscalar);
            p.addParameter('unitDensity', [], @isnumeric);
            
            p.parse(varargin{:});
            
            % set properties
            obj.wave = p.Results.wave(:);
            obj.wave_ = (400:10:700)';
            obj.density = p.Results.density;
            
            if isempty(p.Results.unitDensity)
                obj.unitDensity_ = ieReadSpectra('macularPigment.mat', ...
                    obj.wave_)/0.3521;
            else
                obj.unitDensity = p.Results.unitDensity;
            end
        end
        
        % get methods for dependent variables
        function val = get.unitDensity(obj)
            % interpolate for wavelength samples
            val = interp1(obj.wave_,obj.unitDensity_,obj.wave,'linear',0);
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
        function set.unitDensity(obj, val)
            % interpolate for wavelength samples
            obj.unitDensity_ = interp1(obj.wave,val,obj.wave_,'linear',0);
        end
    end
    
    methods (Static)
        function density = eccDensity(eccDeg, varargin)
            % Compute macular pigment optical density as a function of
            % eccentricity
            %
            %    density = macularDensity(eccentricity)
            %
            % Inputs:
            %   eccDeg - eccentricity in degrees
            %
            % Outputs:
            %   density - macular pigment optical density
            %
            % Notes:
            %   1) Macular pigment density is roughly symmetric and thus we
            %      approximate the 2D position by 1D eccentricity
            %   2) The lorentzian function is fitted from data grabbed from
            %      figure 2(B) in the reference paper. The data is stored
            %      in macularDensity.mat
            %
            % Reference:
            %   Putnam, C. M., & Bland, P. J. (2014). Macular pigment
            %   optical density spatial distribution measured in a subject
            %   with oculocutaneous albinism. Journal of Optometry, 7(4),
            %   241-245.
            
            % Compute density with the lorentz function 
            % 
            % Here, we force the model to be symmetric and have 0 density
            % at infinite eccentricity
            density = 0.35 * 3.6028 ./ (eccDeg.^2 + 3.6028);
        end
    end
end