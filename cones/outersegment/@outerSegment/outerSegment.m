classdef outerSegment < handle
% Parent cl.: Compute cone photocurrent(pA) from isomerization rate(R*)
%
% Syntax:
%   os1 = outerSegment();
%
% Description:
%    This is an abstract class: it is a template but never creates a
%    computable object itself. The computable objects are the subclasses,
%    osLinear and osBioPhys. Those classes store the parameters and methods
%    to compute outer segment current responses in two different ways.
%
%	 * osBioPhys uses the set of difference equations based on physiology
%      from Fred Rieke's lab
% 
%    * osLinear derives a linear temporal filter for small signal
%      deviations from the mean isomerization (R*) rate, again using the
%      Rieke lab model.
%
%    If the noiseFlag property is set to true, the compute methods add
%    noise to the output current. The default noise values are those
%    determined by the Angueyra and Rieke (2013, Nature Neuroscience).
%
% Inputs:
%    None required.
%
% Outputs:
%    The created outersegment object
%
% Optional key/value pairs:
%    None.
%
% References:
%    * <http://isetbio.org/cones/adaptation%20model%20-%20rieke.pdf>
%    * <https://github.com/isetbio/isetbio/wiki/Cone-Adaptation>
%    * <http://www.nature.com/neuro/journal/v16/n11/abs/nn.3534.html>
%
% See Also:
%	 Subclasses: osLinear.m, osBioPhys.m
%
% History:
%    08/xx/15  JRG, NC, DHB Created 8/2015
%    02/12/18  jnm          Formatting

% Public read/write properties
%     properties
%     end

% Public, read-only properties.
properties (SetAccess = private, GetAccess = public)
end

% Protected properties; Methods of the parent class and all of its
% subclasses can set these.
properties
    %noiseFlag - determines whether noise is added  
    noiseFlag;

    %patchSize - spacing between cones (width) in um
    patchSize;

    %eccentricityDegs - eccentricity in degrees
    eccentricityDegs;
end

properties (SetAccess = protected)
    %coneCurrentSignal - output signal in pA
    % coneCurrentSignal;
end

properties (SetObservable, AbortSet)
    %timeStep - sampling interval in sec
    timeStep;
end

properties (Constant)
    %validNoiseFlags - The possible valid noise flags.
    %   Options are 'random', 'frozen', and 'none'
    validNoiseFlags = {'random', 'none', 'frozen'};
end

% Public methods
methods
    function obj = outerSegment(varargin)
        % Initialize defaults for outerSegment parameters
        %
        % Syntax:
        %   obj = outerSegment([varargin]);
        %
        % Description:
        %    Initialize the default values for the noise flag, time step,
        %    and eccentricity degrees parameters. They are 'random' 1e-4,
        %    and 0 degrees respectively.
        %
        % Inputs:
        %    None required.
        %
        % Outputs:
        %    obj - The created outer segment object
        %
        % Optional key/value pairs:
        %    None.
        %
        obj.noiseFlag = 'random';
        % obj.coneCurrentSignal = [];
        obj.timeStep = 1e-4;
        obj.eccentricityDegs = 0;
    end

    function set.noiseFlag(obj, val)
        % Function to assign the noise flag
        %
        % Syntax:
        %   set.noiseFlag(obj, val)
        %
        % Description:
        %    Set the noise flag if the value is a string on the valid noise
        %    flags list.
        %
        % Inputs:
        %    obj - The outer segment object to assign the flag to
        %    val - The string value to assign to noise flag.
        %
        % Outputs:
        %    None.
        %
        % Optional key/value pairs:
        %    None.
        %
        if ischar(val) && ...
                (ismember(lower(val), outerSegment.validNoiseFlags))
            obj.noiseFlag = val;
        else
            s = sprintf('%s ', outerSegment.validNoiseFlags{:});
            error(['''%s'' is an invalid value for os.noiseFlag. ' ...
                'Choose one from: %s '], val, s);
        end
    end

    % see osSet in @osLinear and @osBioPhys for details
    function obj = set(obj, varargin)
        osSet(obj, varargin{:});
    end

    % see osGet in @osLinear and @osBioPhys for details
    function val = get(obj, varargin)
        val = osGet(obj, varargin{:});
    end

end

% Methods that must only be implemented in the subclasses.
% If a subclass does not implement each and every of these methods
% it cannot instantiate objects.
methods (Abstract, Access=public)
    % see osLinearCompute, osBioPhysCompute
    compute(obj, pRate, coneType, varargin);

    % see osLinearPlot, osBioPhysPlot
    % plot(obj, plotType);

end

% Methods may be called by the subclasses, but are otherwise private
methods (Access = protected)    
end

methods (Static)
    resampledPhotocurrents = resample(photocurrents, ...
        originalTimeAxis, resampledTimeAxis);
end

% Methods that are totally private (subclasses cannot call these)
methods (Access = private)
    initialize(obj);
end

end
