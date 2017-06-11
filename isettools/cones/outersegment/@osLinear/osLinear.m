classdef osLinear < outerSegment
    % Linear subclass of the outersegment object to convert isomerizations
    % (R*) to current
    %
    %   os = osLinear();
    %
    % osLinear derives a linear temporal filter for small signal deviations
    % from the mean isomerization (R*) rate.  The linear filters are
    % derived as small signal perturbations of the Rieke lab model.
    %
    %   http://isetbio.org/cones/adaptation%20model%20-%20rieke.pdf
    %   https://github.com/isetbio/isetbio/wiki/Cone-Adaptation
    %
    % osLinear implements isomerizations (R*) to photocurrent (pA) using only
    % cone linear temporal filters.
    %
    % If the noiseFlag property of the osLinear object is set to true, the
    % compute methods add noise to the current output signal.  The default
    % values are those determined by the Angueyra and Rieke (2013, Nature
    % Neuroscience).
    %    <http://www.nature.com/neuro/journal/v16/n11/abs/nn.3534.html>
    % 
    % The current is calculated by convolving separate temporal filters for
    % the L, M and S cones with the isomerization time course.
    %
    % See outerSegment class for more information.
    %
    % JRG/HJ/BW, ISETBIO Team, 2016
    
    % Public, read-only properties.
    properties (GetAccess = public, SetAccess = public)
        % linear filters for LMS cones in columns
        lmsConeFilter;
    end
    
    properties (Access = private)
        % The properties in this part is used to store the state of the
        % outersegment so that we can do incremental computation
        % Will likely be deprecated.
        absHistory = [];  % recent absorption history
    end
    
    % Defined in separate files within the @osLinear directory
    methods (Access=public)
        % Will be deprecated
        % [lmsFilters, meanCurrent] = generateBioPhysFilters(os, meanRate, varargin);
        
        % Calculate linear filters and mean current for a cMosaic
        [lmsFilters, meanCurrent] = linearFilters(os, cMosaic, varargin);
        
    end
    
    % Public methods
    methods
        % Constructor
        function obj = osLinear(varargin)
            % Initialize the parent class
            obj = obj@outerSegment();
        end
        
        % set function, see osLinearSet for details
        function obj = set(obj, varargin)
            osSet(obj, varargin{:});
        end
        
        % get function, see osLinearGet for details
        function val = get(obj, varargin)
            val = osGet(obj, varargin{:});
        end
        
        function val = timeAxis(obj)
            % The temporal samples for the lms filters
            if isempty(obj.lmsConeFilter)
                %warning('No lms impulse response functions computed. Returning empty time axis.');
                val = [];
            else
                % Time axis is the length of the filters multiplied by the
                % time step
                val = ((1:size(obj.lmsConeFilter,1)) - 1) * obj.timeStep;
            end
        end
        
        matchSensor(obj, varargin);
        
        % We need to implement compute because it is part of the abstract
        % outerSegment class.  We should probably rename osCompute to
        % compute (BW).
        function [current, interpFilters, meanCur] = compute(obj, cMosaic, varargin)
         [current, interpFilters, meanCur] = osCompute(obj, cMosaic, varargin{:});
        end
        
        % Declaration for functions in the directory
        plot(obj, pType, varargin);
        
    end

    
    methods (Access = private)
    end
end

