% Class to manage unified and controlled access to the fields of old-style
% and new-style cal files. The reason for this class's existence is to allow
% existing PsychCal functions (e.g., enter name here) to operate on both old-style and new-style cal files
% via a unified parameter access framework.
%
% Example usage:
%
% % Import a cal file
% [cal, ~] = GetCalibrationStructure('Enter calibration filename','ViewSonicProbe',[]);
% 
% % Instantiate @CalStruct object to manage the imported cal
% calStructOBJ = CalStruct(cal, 'verbosity', 1);
% 
% % Get a field (the radiometer's serial no)
% meterSerialNo = calStructOBJ.get('meterSerialNumber');
%
% % Set a field (here, the gammaFormat parameter)
% calStructOBJ.set('gammaFormat', 0);
%
% % Get the modified cal (in old-style format)
% cal = calStructOBJ.cal;
%
% % Pass the modified cal to a PsychCal routine
% ....
% % More examples here
%
% 05/02/2014   npc  Wrote it.
%

classdef CalStruct < handle

    % Read-write properties.
    properties 
        verbosity = 0;
        
    end % Public properties
    
    % Read-only properties
    properties (SetAccess = private) 
        % the old-style cal that is returned to the user
        cal;
        
        % Flag indicating whether the inputCal has new-style format.
        inputCalHasNewStyleFormat;
    end
    
    % invisible - to the user properties
    properties (Access = private)
        
        % The cal struct that we receive during instantiation. 
        % This will be modified via calls to calStruct.set(fieldName, fieldValue);
        inputCal;

        % Dictionary for mapping unified field names
        fieldMap = [];
        
        % Flag indicating whether to report a warning when S_device or
        % S_ambient properties are accessed.
        reportSdev_Samb_warning = true;
        
        % properties holding all the cal fields that can be addressed 
        % by a unified field name. We follow the convention that the
        % unified name is the same as the one used in the old-style cal,
        % whereas the (local) properties follow the new-style cal naming.
        
        % General info
        describe___computerInfo;
        describe___svnInfo;
        describe___matlabInfo;                      % only in new-style cal
        describe___monitor;
        describe___displayDeviceName;
        describe___displayDeviceType;
        describe___calibrationType;                 % only in old-style cal
        describe___promptForName;                   % only in old-style cal
        describe___driver;
        describe___who;
        describe___date;
        describe___executiveScriptName;
        describe___comment;
        
        % Screen configuration
        describe___displaysDescription;
        describe___whichScreen;
        describe___blankOtherScreen;
        describe___whichBlankScreen;
        describe___blankSettings;
        describe___dacsize;
        describe___hz;
        describe___screenSizePixel;
        describe___screenSizeMM;
        
        % HDR-specific
        describe___HDRProjector;                    % only in old-style cal
        describe___yokedmethod;                     % only in old-style cal
        
        % Calibration params
        describe___boxSize;
        describe___boxOffsetX;
        describe___boxOffsetY;
        describe___nAverage;
        describe___nMeas;
        describe___bgColor;
        describe___fgColor;
        describe___displayPrimariesNum;
        describe___primaryBasesNum;
        describe___useBitsPP;
        
        % Gamma-fit params
        describe___gamma;
        describe___gamma___fitType;
        describe___gamma___exponents;
        describe___gamma___contrastThresh;
        describe___gamma___fitBreakThresh;
        describe___gamma___nInputLevels;
        describe___gamma___useweight;
        
        % Radiometer params
        describe___leaveRoomTime;
        describe___meterDistance;
        describe___meterModel;
        describe___meterSerialNumber;               % only in new-style cal
        
        % RawData (gamma measurements)
        rawData___S;
        rawData___gammaInput;
        rawData___gammaTable;
        rawData___gammaCurveSortIndices; 
        rawData___gammaCurveMeasurements; 
        rawData___gammaCurveMeanMeasurements; 
        rawData___gammaCurveMeasurementsInNewFormat;
                 
        % HDR-specific spectra
        rawData___yokedSpectra;
        rawData___yokedSettings;
        
        % RawData (basic linearity measurements)
        basicLinearitySetup___settings;
        rawData___basicLinearityMeasurements1;
        rawData___basicLinearityMeasurements2;
        
        % RawData (background-dependence measurements)
        backgroundDependenceSetup___settings;
        backgroundDependenceSetup___bgSettings;
        rawData___backgroundDependenceMeasurements;
        
        % RawData (ambient measurements)
        rawData___ambientMeasurements;
        
        % ProcessedData
        processedData___monSVs;
        processedData___gammaFormat;
        processedData___gammaInput;
        processedData___gammaTable;
        processedData___S_device;
        processedData___P_device;
        processedData___T_device;
        processedData___S_ambient;
        processedData___P_ambient;
        processedData___T_ambient; 
        
        % RuntimeData
        runtimeData___gammaMode;
        runtimeData___iGammaTable;
        runtimeData___T_sensor;
        runtimeData___S_sensor;
        runtimeData___T_linear;
        runtimeData___S_linear;
        runtimeData___M_device_linear;
        runtimeData___M_linear_device;
        runtimeData___M_linear_device_status;
        runtimeData___M_ambient_linear;
        runtimeData___ambient_linear;
    end
    
    
    % Public methods
    methods
        % Constructor
        function obj = CalStruct(cal, varargin)
            
            parser = inputParser;
			parser.addParameter('verbosity', 0);
            % Execute the parser
			parser.parse(varargin{:});
            
            % Create a standard Matlab structure from the parser results.
			parserResults = parser.Results;
            pNames = fieldnames(parserResults);
            for k = 1:length(pNames)
                obj.(pNames{k}) = parserResults.(pNames{k}); 
            end 
            
            % make a private copy
            obj.inputCal = cal;
            
            % set the field mapping
            obj.setFieldMapping();
            
            % detemine the format (old-style or new-style) of the input cal
            obj.determineInputCalFormat();
            
            % parse the input cal            
            obj.parseInputCal();

        end
        
        % Getter method for cal (always in old-style format)
        function cal = get.cal(obj)
            cal = obj.generateUpdatedCal;
        end
        
        % Getter method for a given unified field name
        fieldValue = get(obj, unifiedFieldName);
        
        % Setter method for a given unified field name
        set(obj, unifiedFieldName, fieldValue);
        
        % Getter methods for select properties over which we want more
        % controlled access
        function value = get.rawData___S(obj)
            % Return S
            value = obj.rawData___S;
            % Check if an 'S_device' field exists and if it does that it's
            % value matches that of S. If not, report an error.
            previous_reportSdev_Samb_warning = obj.reportSdev_Samb_warning;
            obj.reportSdev_Samb_warning = false;
            %fprintf('CalStruct.getS: Will access S_device with no warning\n');
            if (~isempty(obj.processedData___S_device) & ~isnan(obj.processedData___S_device))
                %fprintf('Accessed S_device\n'); pause;
                if (any(value ~= obj.processedData___S_device))
                    fprintf(2,'An ''S_device property exists and does not match ''S''.\n');
                    eval('S_device = obj.processedData___S_device')
                    eval('S = value')
                    error('>>>>>>> Mismatched S and S_device ! <<<<');
                end
            end
            obj.reportSdev_Samb_warning = false;
            %fprintf('CalStruct.getS: Will access S_ambient with no warning\n');
            % Check if an 'S_ambient' field exists and if it does that it's
            % value matches that of S. If not, report an error.
            if (~isempty(obj.processedData___S_ambient) & ~isnan(obj.processedData___S_ambient))
                if (any(value ~= obj.processedData___S_ambient))
                    fprintf(2,'An ''S_ambient property exists and does not match ''S''.\n');
                    eval('S_ambient = obj.processedData___S_ambient')
                    eval('S = value')
                    error('>>>>>>> Mismatched S and S_ambient ! <<<<');
                end
            end
            
            obj.reportSdev_Samb_warning = previous_reportSdev_Samb_warning;
        end
        
        function value = get.processedData___S_device(obj)
            % Return S_device
            value = obj.processedData___S_device;
            if (obj.reportSdev_Samb_warning)
                % warn user that he/she must only use the 'S' variable
                fprintf(2,'Use of the ''S_device'' field is discouraged and will be not allowed in future releases. Use ''S'' instead.\n');         
            end
        end
        
        function value = get.processedData___S_ambient(obj)
            % Return S_ambient
            value = obj.processedData___S_ambient;
            if (obj.reportSdev_Samb_warning)
                % warn user that he/she must only use the 'S' variable
                fprintf(2,'Use of the ''S_ambient'' field is discouraged and will be not allowed in future releases. Use ''S'' instead.\n');         
            end    
        end
        
        % End of getter methods for select properties over which we want more
        % controlled access
        
        % Method to print the field names contained in the FieldMap
        printMappedFieldNames(obj);
        
    end
   
    % Private methods
    methods (Access = private)
        % Method to parse the input cal struct
        parseInputCal(obj);
        
        % Method to determine whether the inputCal has new-style format.
        determineInputCalFormat(obj);
        
        % Method to determine whether the inputCal has the expected basic
        % fields.
        determineInputCalValidity(obj);
        
        % Method to check the validity of the requested unified field name.
        isValid = fieldNameIsValid(obj, unifiedFieldName);
        
        % Method to generate a cal struct with the the old-style format,
        % and which will be returned to routines that rely on that format
        cal = generateUpdatedCal(obj);

        % Method to combine svnInfo and matlabInfo into a single struct
        % as was done in the old-style format
        svn = SVNconversion(obj, propertyValue);
        
        % Method to eliminate the displayDescription of the otherDisplays
        % as was done in the old-style format
        displayDescription = DisplaysDescriptionConversion(obj, propertyValue);
        
        % Method to make old-style meter type. In OOC calibration we store
        % the meter model as a string.
        meterID = MeterTypeConversion(obj, propertyValue);
        
        % Method to generate the old-style monIndex parameter  
        monIndex = MonIndexConversion(obj, propertyValue);
        
        % Method to generate the old-style monSpd parameter  
        monSpd = MonSpdConversion(obj, propertyValue);
        
        % Method to generate the old-style mon parameter  
        mon = MonConversion(obj, propertyValue);
        
        % Method to generate the old-style rawGammaInput parameter  
        rawGammaInput = GammaInputConversion(obj, propertyValue);
        
        % Method to generate the old-style basic linearity spectra parameter
        spectra = BasicLinSpectraConversion(obj, propertyValue);
        
        % Method to generate the old-style background-dependence spectra parameter
        bgSpectra = BGspectraConversion(obj, propertyValue);   
        
    end  % private methods
    
    % Useful static functions
    methods (Static)
        
        % Method to determine if a struct/object contains a field/property with agiven name
        result = isFieldOrProperty(structOrObject, fieldOrPropertyName);
        
    end  % Static functions
end