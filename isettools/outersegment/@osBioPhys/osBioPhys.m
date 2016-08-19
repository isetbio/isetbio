classdef osBioPhys < outerSegment 
% BioPhys subclass of the outersegment object
% 
%       os = osBioPhys();
% 
% Converts isomerizations (R*) to outer segment current (pA). The
% difference equation model by Rieke implements a biophysical
% simulation of the phototransduction cascade. If the noiseFlag
% property of the osLinear object is set to 1, this method will add noise
% to the current output signal.
%
% Reference:
%   http://isetbio.org/cones/adaptation%20model%20-%20rieke.pdf
%   https://github.com/isetbio/isetbio/wiki/Cone-Adaptation
% 
% JRG/HJ/BW, ISETBIO Team, 2016

    properties(Access = public)
        % state;   % biophysics parameter state

        %         % Peripheral parameters
        %             sigma = 22;  % rhodopsin activity decay rate (1/sec) - default 22
        %             phi = 22;     % phosphodiesterase activity decay rate (1/sec) - default 22
        %             eta = 2000;	  % phosphodiesterase activation rate constant (1/sec) - default 2000
        %             gdark = 20.5; % concentration of cGMP in darkness - default 20.5
        %             k = 0.02;     % constant relating cGMP to current - default 0.02
        %             h = 3;       % cooperativity for cGMP->current - default 3
        %             cdark = 1;  % dark calcium concentration - default 1
        %             beta = 9;	  % rate constant for calcium removal in 1/sec - default 9
        %             betaSlow = 0.4; % rate constant for slow calcium modulation of channels - default 0.4
        %             n = 4;  	  % cooperativity for cyclase, hill coef - default 4
        %             kGc = 0.5;   % hill affinity for cyclase - default 0.5
        %             OpsinGain = 10; % so stimulus can be in R*/sec (this is rate of increase in opsin activity per R*/sec) - default 10
        

        % Foveal parameters
        sigma = 10;       % rhodopsin activity decay rate (1/sec) - default 22
        phi   = 22;       % phosphodiesterase activity decay rate (1/sec) - default 22
        eta   = 700;      % phosphodiesterase activation rate constant (1/sec) - default 2000
        gdark = 20.5;     % concentration of cGMP in darkness - default 20.5
        k     = 0.02;     % constant relating cGMP to current - default 0.02
        h     = 3;        % cooperativity for cGMP->current - default 3
        cdark = 1;        % dark calcium concentration - default 1
        beta  = 5;        % rate constant for calcium removal in 1/sec - default 9
        betaSlow = 0.4;   % rate constant for slow calcium modulation of channels - default 0.4
        n     = 4;        % cooperativity for cyclase, hill coef - default 4
        kGc   = 0.5;      % hill affinity for cyclase - default 0.5
        OpsinGain = 12;   % so stimulus can be in R*/sec (this is rate of increase in opsin activity per R*/sec) - default 10
        
    end
    

    
    methods
        function obj = set(obj, varargin)
            % set function, see osBioPhysSet for details
            osSet(obj, varargin{:});
        end
        
        
        function val = get(obj, varargin)
            % get function, see osBioPhysGet for details
            val = osGet(obj, varargin{:});
        end
    end
    
    methods (Access=public)        
        function obj = compute(obj, sensor, varargin)
            % see osCompute for details
            obj = osCompute(obj, sensor, varargin{:});
        end
        
        % Derive some parameters - steady state constraints among parameters
        function q = get.q
            q    = 2 * beta * cdark / (k * gdark^h);
        end
        
        function smax = get.smax
            smax = eta/phi * gdark * (1 + (cdark / kGc)^n);
        end
    
        function plot(obj, sensor, varargin)
            % see osPlot for details
            osPlot(obj, sensor, varargin{:});
        end
    end
end
