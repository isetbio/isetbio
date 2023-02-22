classdef RGCRF < handle

    % A class to handle the spatiotemporal RF dynamics of RGCs
    %
    % Syntax:
    %   
    % Description
    %

    % Public properties
    properties  (GetAccess=public, SetAccess=public)
        
    end % Public properties

    % Read-only properties
    properties (GetAccess=public, SetAccess=private)
        % The class type
        classType;

        % The spatial RF
        spatialRF

        % The temporal RF
        temporalRF

        % The temporal frequency support
        tfSupport = 0:0.5:500;

        % Container with temporal kernel params for various RGC types
        temporalKernelParams;
    end

    % Constant properties
    properties (Constant)
        validClassTypes = {
            'midget' ...
            };

        % Temporal parameters for the ON midget cell example from Fig 6 of
        % Benardete & Kaplan '97: "The receptive field of the primate P retinal ganglion cell, I: Linear dynamics"
        midgetON_fig6Example_Pstruct = struct(...
            'center', struct(...)
                'gainImpulsesPerSecondPerUnitContrast', 184.2, ...
                'gSubtractive', 0.69, ...
                'tauHighPassMsec', 18.61,...
                'lowPassStagesNum', 38, ...
                'tauLowPassMsec', 1.23, ...
                'initialDelayMsec', 4.0 ...
                ), ...
            'surround', struct(...)
                'gainImpulsesPerSecondPerUnitContrast', -125.33, ...
                'gSubtractive', 0.56, ...
                'tauHighPassMsec', 33.28,...
                'lowPassStagesNum', 124, ...
                'tauLowPassMsec', 0.42, ...
                'initialDelayMsec', 4.0 ...
                ) ...
            );

        midget_median_Pstruct = struct(...
            'center', struct(...)
                'gainImpulsesPerSecondPerUnitContrast', 67.59, ...
                'gSubtractive', 0.69, ...
                'tauHighPassMsec', 29.36,...
                'lowPassStagesNum', 38, ...
                'tauLowPassMsec', 48.15/38, ...
                'initialDelayMsec', 3.5 ...
                ), ...
            'surround', struct(...)
                'gainImpulsesPerSecondPerUnitContrast', 49.98, ...
                'gSubtractive', 0.48, ...
                'tauHighPassMsec', 18.62,...
                'lowPassStagesNum', 111, ...
                'tauLowPassMsec', 55.02/111, ...
                'initialDelayMsec', 3.5 ...
                ) ...
            );


        % Temporal parameters for the OFF midget cell example from Fig 6 of
        % Benardete & Kaplan '97: "The receptive field of the primate P retinal ganglion cell, I: Linear dynamics"
        midgetOFF_fig6Example_Pstruct = struct(...
            'center', struct(...)
                'gainImpulsesPerSecondPerUnitContrast', -114.12, ...
                'gSubtractive', 0.82, ...
                'tauHighPassMsec', 24.9,...
                'lowPassStagesNum', 25, ...
                'tauLowPassMsec', 2.12, ...
                'initialDelayMsec', 3.5 ...
                ), ...
            'surround', struct(...)
                'gainImpulsesPerSecondPerUnitContrast', 74.57, ...
                'gSubtractive', 0.72, ...
                'tauHighPassMsec', 49.81,...
                'lowPassStagesNum', 83, ...
                'tauLowPassMsec', 0.76, ...
                'initialDelayMsec', 3.5 ...
                ) ...
            );
         
    end

    % Public methods
    methods

        % Constructor
        function obj = RGCRF(classType, varargin)
            % Validate the class type
            assert(ismember(classType, RGCRF.validClassTypes), '''%s'' is not a valid RGCRF class type.', classType);
            obj.classType = classType;

            % Parse optional input
            p = inputParser;
            p.addParameter('visualizeDynamics', false, @islogical);
            p.parse(varargin{:});

            visualizeDynamics = p.Results.visualizeDynamics;

            % Select pStructs for this class type
            switch (classType)
                case 'midget'
                    ids = {'midgetON_fig6Example', 'midgetOFF_fig6Example', 'midgetMedian'};
                    pStructs = {RGCRF.midgetON_fig6Example_Pstruct, RGCRF.midgetOFF_fig6Example_Pstruct, RGCRF.midget_median_Pstruct};
            end

            % Populate container with pStructs
            obj.temporalKernelParams = containers.Map(ids, pStructs);
            
            cellID = 'midgetON_fig6Example';
            obj.generateMidgetRF(cellID, visualizeDynamics);

        end % Constructor

    end % Public methods


    % Private methods
    methods (Access=private)
        % Method to generate the RF dynamics of a midget RGC
        generateMidgetRF(obj, cellID, visualizeDynamics);

        % Equation (2) of Benardete & Kaplan '97:
        % "The receptive field of the primate P retinal ganglion cell, I: Linear dynamics"
        temporalTransferFunction = generateTemporalTrasferFunction(obj, p);
    end % Private methods


    % Static methods
    methods (Static)

        % Generate impulse response from one-sided temporal transfer function
        [impulseResponse, temporalSupport] = impulseResponseFromOneSidedTransferFunction(oneSidedTransferFunction, tfSupport)

    end % Static methods

end