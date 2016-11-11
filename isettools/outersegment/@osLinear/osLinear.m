classdef osLinear < outerSegment
% Linear subclass of the outersegment object
%
%   os = osLinear();
%
% Implements isomerizations (R*) to photocurrent (pA) using only cone
% linear temporal filters. The default values are those determined by
% the Angueyra and Rieke (2013, Nature Neuroscience).
%
% The current is calculated by convolving separate temporal filters
% for the L, M and S cones with the isomerization time course.
%
% JRG/HJ/BW, ISETBIO Team, 2016
    
    % Public, read-only properties.
    properties (GetAccess = public, SetAccess = private)
        % linear filters for LMS cones in columns
        lmsConeFilter;
    end
    
    properties (Access = private)
        % The properties in this part is used to store the state of the
        % outersegment so that we can do incremental computation
        
        absHistory = [];  % recent absorption history 
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
        
        matchSensor(obj, varargin);
        
        function obj = compute(obj, pRate, coneType, varargin)
            % see osCompute for details
            obj = osCompute(obj, pRate, coneType, varargin{:});
        end
        
        function plot(obj, cMosaic)
            % see osPlot for details
            osPlot(obj, cMosaic);
        end
    end
    
    methods (Access = private)
        function newIRFs = generateLinearFilters(obj, meanRate, varargin)
            % Generates the linear temporal filters for the L, M and S
            % cones based on physiological data from Angueyra and Rieke
            % (Nat. Neuro., 2013).
            %
            %   newIRFs = generateLinearFilters([coneMosaic])
            %
            % Inputs: 
            %   obj      - osLinear class object
            %   meanRate - mean photon isomerization rate
            %
            % Outputs:
            %  newIRFs  - nx3 matrix with each row the temporal impulse
            %             response for the L, M and S cones, respectively.
            %
            % See also:
            %   osFitlerConesLinear
            %
            % JRG/HJ/BW, ISETBIO TEAM, 2016
            
            % Handle input
            p = inputParser;
            p.addRequired('obj', @(x) isa(x, 'osLinear'));
            p.addParameter('osType',0,@islogical);
            p.addRequired('meanRate', @isscalar);
            
            p.parse(obj, meanRate, varargin{:});          
                        
            osType = p.Results.osType; % peripheral (0) or foveal (1)
            
            dt = obj.timeStep;
            timeAxis = 0 : dt : 0.3;
            
            % parameters of cone temporal filters measured in Angueyra and
            % Rieke (Nat. Neuro., 2013).
            scFact = 1;      % scale factor
            tauR = 0.0216;   % rising phase time constant
            tauD = 0.0299;   % damping time constant
            tauP = 0.5311;   % period
            phi = 34.1814;   % phase
            
            % generate filter
            filter = scFact * (timeAxis/tauR).^3 ...
                ./ (1 + (timeAxis/tauR).^3) .* exp(-(timeAxis./tauD)) ...
                .* cos(2*pi*timeAxis/tauP+ deg2rad(phi));
            
            % Adjust gain of the cone impulse response for L, M, S cones
            %
            % Weber-Fechner relation:
            %  gain/gain_dark = ...
            %    1/(1+(bgIntensity/halfDesensitizationIntensitity))
            %            
            Io = 2250; % half-desensitizing background in R*/cone/sec
            Ib = [7131 6017 1973]; % R* per sec due to background adapting
            
            % Set gain parameters according to Angueyra & Rieke:
            Ib = Ib * meanRate / max(Ib);
            gain_dark = 0.22;              % from Juan's paper
            
            % scale IRF to reflect amplitude at chosen background
            % using Weber adaptation equation
            gain = gain_dark ./ (1 + Ib./Io);
            
            % The L, M and S IRFs are created by scaling Filter with the
            % scale factors
            newIRFs = filter(:) * gain * dt ./ max(filter);
            obj.lmsConeFilter = newIRFs;
        end
        
        function impulseResponseLMS = generateBioPhysFilters(os, varargin)
            % Generates the impulse responses for the LMS cones with the appropriate
            % magnitude using the biophysical model for the outer segment.
            %
            % Manually sets an impulse stimulus in the absorptions field of a new cone
            % mosaic object, generates the impulse response for the LMS cones and
            % returns them, in order to be used in a linear computation with a
            % near-constant background.
            %
            % See also: v_osBioPhys, t_coneMosaicFoveal, t_osLinearize
            %
            % 11/2016 JRG (c) isetbio team
            
            % parse input parameters
            p = inputParser; p.KeepUnmatched = true;
            p.addRequired('os', @(x) isa(x, 'outerSegment'));
            p.addParameter('meanRate',@isnumeric);
            p.addParameter('osType',false,@islogical);
            p.addParameter('coneType',@isnumeric);
            p.parse(os, varargin{:})
            
            os = p.Results.os;
            osType = p.Results.osType;
            meanIsoArray = p.Results.meanRate;
            coneType = p.Results.coneType;
            
            % Generate impulse responses for L, M or S cone
            % A new cone mosaic is generated for each type of cone, and the
            % impulse response
            
            timeStep = os.timeStep;        % time step
            nSamples = round(0.8/timeStep) + 1;       % 0.3 sec
            
            % Manually set flash rate to just one at 5000 R*/sec
            flashIntens = 10;         
            
            backgroundMean = 1000;   %[10 100 500 1000 2000 5000 10000];
            
            flashTimeStep = 5000;    % Warm up period, then flash
            
%             vcNewGraphWin; hold on           
            
            for meanInd = 1:length(meanIsoArray)
                
                % Build the stimulus
                stimulus = ones(1,nSamples)*backgroundMean;
                
                % Set mean background for start of current output
                bgIsomerizations = backgroundMean;
                
                % flash intensity in R*/cone/sec (maintained for 1 bin only)
                meanIntens  = meanIsoArray(meanInd);
                
                % Create stimulus.
                stimulus = meanIntens*ones(nSamples, 1);
                stimulusFlat = reshape(stimulus, [1 1 nSamples]);
                
                % stimulus(1) = flashIntens*timeStep;
                stimulus(flashTimeStep) = flashIntens;
                stimulus = reshape(stimulus, [1 1 nSamples]);
                
                % Generate the cone mosaic for the impulse stimulus
                osCM = osBioPhys('osType',osType);            
                osCM.set('noise flag',0);
                cm = coneMosaic('os',osCM,'pattern', 2); % a single cone
                cm.integrationTime = timeStep;
                cm.os.timeStep = timeStep;
                
                cm.absorptions  = stimulus;
                
                % Compute outer segment currents.
                cm.computeCurrent('bgR',bgIsomerizations);
                currentScaled = (cm.current);% - cm.current(flashTimeStep);
                
                % Generate the cone mosaic for the flat stimulus
                osCMflat = osBioPhys('osType',osType);            
                osCMflat.set('noise flag',0);
                cmFlat = coneMosaic('os',osCMflat ,'pattern', 2); % a single cone
                cmFlat.integrationTime = timeStep;
                cmFlat.os.timeStep = timeStep;
                
                cmFlat.absorptions  = stimulusFlat;
                
                % Compute outer segment currents.
                cmFlat.computeCurrent('bgR',bgIsomerizations);
                currentScaledFlat = (cmFlat.current);% - cm.current(flashTimeStep);
                
                % Plot the impulse responses for comparison.
                
%                 subplot(1,2,1); hold on;
%                 tme = (1:nSamples)*timeStep;
%                 plot(tme,squeeze(currentScaled),'LineWidth',3);
%                 grid on; % legend('Peripheral','Foveal');
%                 % axis([0 0.2 min((current2Scaled./max(current2Scaled(:)))) 1]);
%                 xlabel('Time (sec)','FontSize',14);
%                 ylabel('pA','FontSize',14);
%                 set(gca,'fontsize',14);
%                 title(sprintf('Peripheral (fast)\nImpulse response as a function of mean R*/sec'));
%                 
%                 subplot(1,2,2); hold on;
%                 tme = (1:nSamples)*timeStep;
%                 plot(tme-flashTimeStep*timeStep,squeeze(currentScaled)- squeeze(currentScaledFlat),'LineWidth',3);
%                 grid on; % legend('Peripheral','Foveal');
%                 % axis([0 0.2 min((current2Scaled./max(current2Scaled(:)))) 1]);
%                 axis([0 nSamples*timeStep-flashTimeStep*timeStep -.2 2]);
%                 xlabel('Time (sec)','FontSize',14);
%                 ylabel('pA','FontSize',14);
%                 set(gca,'fontsize',14);
%                 % title('Impulse response in the dark
%                 
%                 title(sprintf('Peripheral (fast)\nIR as a function of mean R*/sec, rezeroed'));
                
                numberIsomerizations = flashIntens;
                impulseResponseLMS(:,meanInd) = (squeeze(currentScaled(flashTimeStep:end))- squeeze(currentScaledFlat(flashTimeStep:end)))./numberIsomerizations;
            end
            
        end
    end
end
