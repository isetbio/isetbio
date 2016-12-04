classdef outerSegment < handle
    % outerSegment parent class for computing cone photocurrent (pA) from
    % isomerization rate (R*)
    %
    % This is an abstract class: it is a template but never creates a
    % computable object itself. The computable objects are the subclasses,
    % osLinear and osBioPhys.  Those classes store the parameters and
    % methods to compute outer segment current responses in two different
    % ways.
    %
    %  * osBioPhys uses the set of difference equations based on physiology
    % from Fred Rieke's lab
    % 
    %  * osLinear derives a linear temporal filter for small signal
    % deviations from the mean isomerization (R*) rate, again using the
    % Rieke lab model.
    %
    %   <http://isetbio.org/cones/adaptation%20model%20-%20rieke.pdf>
    %   <https://github.com/isetbio/isetbio/wiki/Cone-Adaptation>
    %
    % If the noiseFlag property is set to true, the compute methods add
    % noise to the output current.  The default noise values are those
    % determined by the Angueyra and Rieke (2013, Nature Neuroscience).
    %
    %   <http://www.nature.com/neuro/journal/v16/n11/abs/nn.3534.html>
    %
    % See the subclasses:
    %       osLinear.m, osBioPhys.m
    %
    % JRG, NC, DHB, 8/2015
    
    %     % Public read/write properties
    %     properties
    %     end
    
    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
    end
    
    % Protected properties; Methods of the parent class and all of its
    % subclasses can set these.
    properties
        noiseFlag;            % determines whether noise is added
        patchSize;            % spacing between cones (width) in um
    end
    
    properties (SetAccess = protected)
        % coneCurrentSignal;    % output signal in pA
    end
    
    properties (SetObservable, AbortSet)
        timeStep;             % sampling interval in sec
    end
    
    % Public methods
    methods
        function obj = outerSegment(varargin)
            obj.noiseFlag = 'random';
            % obj.coneCurrentSignal = [];
            obj.timeStep = 1e-4;
        end
        
        % see osSet in @osLinear and @osBioPhys for details
        function obj = set(obj, varargin)
            osSet(obj, varargin{:});
        end
        
        % see osGet in @osLinear and @osBioPhys for details
        function val = get(obj, varargin)
            val = osGet(obj, varargin{:});
        end
        
        function val = timeAxis(obj)
            % The temporal samples for the lms filters
            if isempty(obj.lmsConeFilter)
                error('No filters computed');
            else
                val = ((1:size(obj.lmsConeFilter,1)) - 1) * obj.timeStep;
            end
            
        end
        
    end
    
    % Methods that must only be implemented in the subclasses.
    % If a subclass does not implement each and every of these methods
    % it cannot instantiate objects.
    methods (Abstract, Access=public)
        % see osLinearCompute, osBioPhysCompute
        compute(obj, pRate, coneType, varargin);
        % see osLinearPlot, osBioPhysPlot
        plot(obj, plotType);
        
    end
    
    % Methods may be called by the subclasses, but are otherwise private
    methods (Access = protected)
        
    end
    
    methods (Static)
        resampledPhotocurrents = resample(photocurrents, originalTimeAxis, resampledTimeAxis);
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
        initialize(obj);
    end
    
end
