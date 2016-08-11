classdef Cone < hiddenHandle
    % Class for single cone
    %
    %   cone = Cones();
    %
    % This class contains properties for single cone cell. For the full
    % cone mosaic, see coneMosaic class
    %
    % In this class, we assume cones are squares with same height and width
    %
    % HJ, ISETBIO Team, 2016
    
    properties  % public properties
        opticalDensity;  % photopigment optical densities for L,M,S
        peakEfficiency;  % peak absorptance efficiency
        diameter;        % cone aperture diameter in meters
    end
    
    properties (SetObservable, AbortSet)
        wave;            % wavelength samples
    end
    
    properties (Dependent)
        absorbance;         % spectral absorbance of the cones
        absorptance;        % cone absorptance without ocular media
        quantaFundamentals; % normalized cone absorptance
        energyFundamentals; % normalized cone absorption in energy units
        area;               % diameter^2
    end
    
    properties(Access = private)  % private properties
        absorbanceF_;  % cell array of LMS interpolant for absorbance data
    end
    
    methods  % public methods
        % constructor
        function obj = Cone(varargin)
            % parse input
            p = inputParser;
            p.addParameter('wave', 400:10:700, @isnumeric);
            p.addParameter('opticalDensity', [0.5 0.5 0.4], @isnumeric);
            p.addParameter('absorbance', [], @isnumeric);
            p.addParameter('peakEfficiency', [2 2 2]/3, @isnumeric);
            p.addParameter('diameter', 2e-6, @isnumeric);
            
            p.parse(varargin{:});
            
            % set object properties
            obj.wave = p.Results.wave(:);
            obj.opticalDensity = p.Results.opticalDensity(:);
            obj.peakEfficiency = p.Results.peakEfficiency(:);
            obj.diameter = p.Results.diameter;
            
            if isempty(p.Results.absorbance)
                [logAbs, w] = ieReadSpectra('coneAbsorbance');
                obj.absorbanceF_ = cell(1, 3);
                for ii = 1 : 3
                    % build interpolant for LMS
                    obj.absorbanceF_{ii} = griddedInterpolant(w, ...
                        10.^logAbs(:, ii));
                end 
            else
                obj.absorbance = p.Results.absorbance;
            end
        end
        
        function str = description(obj, varargin)
            % generate description string for this object
            str = sprintf('\tDiameter:\t %.2f um\n', obj.diameter*1e6);
            str = [str sprintf('\tOptical density: ') ...
                sprintf('[%.2f, %.2f, %.2f]\n', obj.opticalDensity(1), ...
                obj.opticalDensity(2), obj.opticalDensity(3))];
            str = [str sprintf('\tPeak efficiency: ') ...
                sprintf('[%.2f, %.2f, %.2f]\n', obj.peakEfficiency(1), ...
                obj.peakEfficiency(2), obj.peakEfficiency(3))];
        end
        
        % get method for dependent variable
        function val = get.absorbance(obj) % inerpolate for absorbance
            val = cell2mat(cellfun(@(x) x(obj.wave), obj.absorbanceF_, ...
                'UniformOutput', false));
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
            val = obj.diameter^2;
        end
        
        % set method for dependent variable
        function set.absorbance(obj, val)
            for ii = 1 : 3
                obj.absorbanceF_{ii} = griddedInterpolant(obj.wave, ...
                    val(:, ii));
            end
        end
        
        function set.area(obj, val)
            obj.diameter = sqrt(val);
        end
    end
end