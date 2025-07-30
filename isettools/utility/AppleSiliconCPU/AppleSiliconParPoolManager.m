classdef AppleSiliconParPoolManager < handle

    % Create an AppleSiliconParPoolManager
    %
    % Syntax:
    %   (Default ParPool manager: Configure a parpool that gives at least 4 GB RAM to each core, and do not print anything on the console)
    %   M3ParPoolManager = AppleSiliconParPoolManager('runSilent', true);
    %   
    %   (Default ParPool manager: Configure a parpool that gives at least 4 GB RAM to each core)
    %   M3ParPoolManager = AppleSiliconParPoolManager('default');
    %  
    %   (Conservative ParPool manager: Configure a parpool that gives at least 8 GB RAM to each core)
    %   M3ParPoolManager = AppleSiliconParPoolManager('consevative');
    %  
    %   (Extreme ParPool manager: Configure a parpool with the max # of cores, independent of system RAM)
    %   M3ParPoolManager = AppleSiliconParPoolManager('extreme');
    % 
    %   (Custom ParPool manager: Configure a parpool with the desired # of cores, independent of system RAM)
    %   nCores = some number;
    %   M3ParPoolManager = AppleSiliconParPoolManager(nCores);
    %  
    %   Restore previous parpool size:
    %   M3ParPoolManager.restoreLastParpoolSize();
    %   
    %   Info on current Manager: Returns a struct with all properties
    %   M3ParPoolManager.info
    %
    % Description:
    %    An object to optimize parpool creation on Apple Silicon machines
    %
    % History:
    %    December 2023  NPC  Wrote it

    
    % Public properties
    properties
        runSilent = false;
    end

    % Read-only properties
    properties (GetAccess=public, SetAccess=private)
        currentParpoolSize = [];
        lastParpoolSize = [];
        parpoolConfig = [];
        currentRAMperCore;
        usablePhysicalMemoryGB;
        physicalMemoryGB;
        coresNum;
        architecture;
    end

    % Private properties
    properties (GetAccess=private, SetAccess=private)
        minRAMperCore;
    end

    % Public methods
    methods
        
        % Constructor
        function obj = AppleSiliconParPoolManager(parpoolConfig, varargin)

            if (nargin == 0)
                obj.parpoolConfig = 'default';
                obj.setParpoolSize();
                return;
            end

            parpoolConfigErrorMessage = 'AppleSiliconParPool() can be initialized either with no arguments, or a single argument which must be either a scalar (>=1), ''default'', ''conservative'', ''extreme''.'; 

            if ~isempty(varargin) && ~isstruct(varargin{1})
                varargin = ieParamFormat(varargin);
            end

            % Parse input
            p = inputParser;
            p.addRequired('parpoolconfig', @(x)parpoolConfigValidationFunction(x, parpoolConfigErrorMessage ));
            p.addParameter('runsilent', false, @islogical);
            p.parse(parpoolConfig, varargin{:});
            obj.parpoolConfig = p.Results.parpoolconfig;
            obj.runSilent = p.Results.runsilent;

            obj.setParpoolSize();

            % Validation function for parpoolConfig
            function s = parpoolConfigValidationFunction(x, parpoolConfigErrorMessage)
                c(1) = ~isempty(x);
                c(2) = ischar(x)&&(~ismember(x, {'default', 'conservative', 'extreme', 'half max'}));
                c(3) = isscalar(x)&&(~(x >=1));
                if (all(c))
                    error(parpoolConfigErrorMessage);
                else
                   s = true;
                end
            end
        end

        function restoreLastParpoolSize(obj)
            if (~isempty(obj.lastParpoolSize))
                % If no pool, do not create new one.
                poolobj = gcp('nocreate'); 
                if (~isempty(poolobj))
                    delete(poolobj);
                end
                if (~obj.runSilent)
                    fprintf('Staring new parpool with the previous # of workers (%d)\n', obj.lastParpoolSize);
                end
                parpool(obj.lastParpoolSize);
            else
                if (~obj.runSilent)
                    fprintf('Did not detect a previously running parpool. Not restoring a previous parpool size.\n');
                end
            end
        end

        function s = info(obj)
            props = properties(obj);
            for iProp = 1:numel(props)
                fieldname = props{iProp};
                s.(fieldname) = obj.(fieldname);
            end
        end

    end % Public methods

    % Private methods
    methods (Access=private)

        function setParpoolSize(obj)
            determinedSystemResourcesSuccessfully = obj.cpuInfo();
            
            if (~determinedSystemResourcesSuccessfully)
                fprintf('Could not determine system resources. No change in parpool.\n');
            end

            if (strcmp(obj.architecture, 'arm'))
                if (isscalar(obj.parpoolConfig))
                    numWorkers = obj.parpoolConfig;
                    obj.minRAMperCore = [];
                    obj.restartParpool(numWorkers);

                elseif (ischar(obj.parpoolConfig))
                    switch (obj.parpoolConfig)
                        case 'default'
                            % Mini GB per core
                            obj.minRAMperCore = 4.0;
                        case 'conservative'
                            obj.minRAMperCore = 8.0;
                        case {'extreme', 'half max', '3/4'}
                            obj.minRAMperCore = obj.usablePhysicalMemoryGB/obj.coresNum;

                        otherwise
                            fprintf('\n\n----> Valid configs are either a scalar between %d and %d, or {''default'', ''conservative'', ''extreme''}. <----\n\n', 1, obj.coresNum);
                            if (ischar(obj.parpoolConfig))
                                 error('''%s'' is an invalid parpoolConfig!', obj.parpoolConfig);
                            else
                                error('%g is an invalid parpoolConfig!', obj.parpoolConfig);
                            end

                    end % switch

                    numWorkers = min([obj.coresNum , max([1 floor(obj.usablePhysicalMemoryGB / obj.minRAMperCore)])]);

                    if (strcmp(obj.parpoolConfig, 'half max'))
                        numWorkers = ceil(0.5*numWorkers);
                    end

                    
                    if (strcmp(obj.parpoolConfig, '3/4'))
                        numWorkers = floor(0.75*numWorkers);
                    end

                    obj.restartParpool(numWorkers);
                else
                    obj.parpoolConfig
                    error('Unknown parpoolConfig!')
                end
            else
                if (~obj.runSilent)
                    fprintf('Not an Apple Silicon machine. No change in parpool.\n');
                end
            end
        end % setParpoolSize

        function determinedSystemResourcesSuccessfully = cpuInfo(obj)

            determinedSystemResourcesSuccessfully = true;

            % Determine system architecture
            [~, systemArchitecture] = system('uname -p');
            if (isempty(systemArchitecture))
                determinedSystemResourcesSuccessfully = false;
                return;
            end
            obj.architecture = systemArchitecture(1:end-1);
    

            % Determine # of cores
            [~, numCPUs] = system('sysctl -n hw.ncpu');
            if (isempty(numCPUs))
                determinedSystemResourcesSuccessfully = false;
                return;
            end
            obj.coresNum = str2double(numCPUs(1:end-1));
            
            [~, physicalMemory] = system('sysctl -a -h | grep memsize:');
            if (isempty(physicalMemory))
                obj.physicalMemoryGB = [];
            else
                physicalMemory = physicalMemory(2:end-1);
                physicalMemory = strrep(physicalMemory, 'w.memsize:', '');
                obj.physicalMemoryGB = str2double(physicalMemory)/1e9;
            end

            [~, usablePhysicalMemory] = system('sysctl -a -h | grep memsize_usable:');
            if (isempty(usablePhysicalMemory))
                obj.usablePhysicalMemoryGB = [];
            else
                usablePhysicalMemory = usablePhysicalMemory(2:end-1);
                usablePhysicalMemory = strrep(usablePhysicalMemory, 'w.memsize_usable:', '');
                obj.usablePhysicalMemoryGB = str2double(usablePhysicalMemory)/1e9;
            end

            if (isempty(obj.usablePhysicalMemoryGB))
                obj.usablePhysicalMemoryGB = obj.physicalMemoryGB;
            end

            if (isempty(obj.physicalMemoryGB))
                obj.physicalMemoryGB = obj.usablePhysicalMemoryGB;
            end

            if (isempty(obj.usablePhysicalMemoryGB)) && (isempty(obj.physicalMemoryGB))
                determinedSystemResourcesSuccessfully = false;
            end

        end

        function restartParpool(obj, desiredNumWorkers)

            if (desiredNumWorkers < 1) || (desiredNumWorkers > obj.coresNum)
                error('Number of workers must be between 1 and %d. %d is not valid.\n', obj.coresNum, desiredNumWorkers);
            end

            % If no pool, do not create new one.
            poolobj = gcp('nocreate');
            if (~isempty(poolobj)) && ((poolobj.NumWorkers == desiredNumWorkers))
                obj.currentParpoolSize = desiredNumWorkers;
                obj.currentRAMperCore = obj.usablePhysicalMemoryGB / obj.currentParpoolSize;
                if (~obj.runSilent)
                    fprintf('Previously established parpool had the desired num of workers. Keeping it as is.\n');
                end
                return;
            end

            if (~isempty(poolobj)) && ((poolobj.NumWorkers ~= desiredNumWorkers))
                if (~obj.runSilent)
                    fprintf('Closing old parpool which had %d workers\n', poolobj.NumWorkers);
                end
                obj.lastParpoolSize = poolobj.NumWorkers;
                delete(poolobj);
            else
                obj.lastParpoolSize = [];
            end

            if (isempty(obj.lastParpoolSize)) || (obj.lastParpoolSize ~= desiredNumWorkers)
                if (~obj.runSilent)
                    fprintf('Starting new parpool with %d workers\n', desiredNumWorkers);
                end
                parpool(desiredNumWorkers);
                obj.currentParpoolSize = desiredNumWorkers;
                obj.currentRAMperCore = obj.usablePhysicalMemoryGB / obj.currentParpoolSize;
            end
        end

    end % Private methods


    % Static methods
    methods (Static)
    end
end