classdef AppleSiliconParPoolManager < handle

    % Create an AppleSiliconParPoolManager
    %
    % Syntax:
    %   (Default ParPool manager: Configure a parpool that gives at least 4 GB RAM to each core)
    %   M3ParPoolManager = AppleSiliconParPoolManager;
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

    % Description:
    %    An object to optimize parpool creation on Apple Silicon machines
    %
    % History:
    %    December 2023  NPC  Wrote it

    
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
        function obj = AppleSiliconParPoolManager(parpoolConfig)

            if (nargin == 0)
                obj.parpoolConfig = 'default';
                obj.setParpoolSize();
                return;
            end

            parpoolConfigErrorMessage = 'AppleSiliconParPool() can be initialized either with no arguments, or a single argument which must be either a scalar (>=1), ''default'', ''conservative'', ''extreme''.'; 

            % Parse input
            p = inputParser;
            p.addRequired('parpoolConfig', @(x)parpoolConfigValidationFunction(x, parpoolConfigErrorMessage ));
            p.parse(parpoolConfig);
            obj.parpoolConfig = p.Results.parpoolConfig;

            obj.setParpoolSize();

            % Validation function for parpoolConfig
            function s = parpoolConfigValidationFunction(x, parpoolConfigErrorMessage)
                c(1) = ~isempty(x);
                c(2) = ischar(x)&&(~ismember(x, {'default', 'conservative', 'extreme'}));
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
                fprintf('Staring new parpool with the previous # of workers (%d)\n', obj.lastParpoolSize);
                parpool(obj.lastParpoolSize);
            else
                fprintf('Did not detect a previously running parpool. Not restoring a previous parpool size.\n')
            end
        end

    end % Public methods

    % Private methods
    methods (Access=private)

        function setParpoolSize(obj)
            obj.cpuInfo();
            
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
                        case 'extreme'
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
                    obj.restartParpool(numWorkers);
                else
                    obj.parpoolConfig
                    error('Unknown parpoolConfig!')
                end
            else
                fprintf('Not an Apple Silicon machine. No change in parpool.\n')
            end
        end % setParpoolSize

        function cpuInfo(obj)
            [~, architecture] = system('uname -p');
            architecture = architecture(1:end-1);
            obj.architecture = architecture;
    
            [~, numCPUs] = system('sysctl -n hw.ncpu');
            obj.coresNum = str2double(numCPUs(1:end-1));
            
                [~, physicalMemory] = system('sysctl -a -h | grep memsize:');
                physicalMemory = physicalMemory(2:end-1);
                physicalMemory = strrep(physicalMemory, 'w.memsize:', '');
                obj.physicalMemoryGB = str2num(physicalMemory)/1e9;
    
                [~, usablePhysicalMemory] = system('sysctl -a -h | grep memsize_usable:');
                usablePhysicalMemory = usablePhysicalMemory(2:end-1);
                usablePhysicalMemory = strrep(usablePhysicalMemory, 'w.memsize_usable:', '');
                obj.usablePhysicalMemoryGB = str2num(usablePhysicalMemory)/1e9;
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
                fprintf('Previously established parpool had the desired num of workers. Keeping it as is.\n');
                return;
            end

            if (~isempty(poolobj)) && ((poolobj.NumWorkers ~= desiredNumWorkers))
                fprintf('Closing old parpool which had %d workers\n', poolobj.NumWorkers);
                obj.lastParpoolSize = poolobj.NumWorkers;
                delete(poolobj);
            else
                obj.lastParpoolSize = [];
            end

            if (isempty(obj.lastParpoolSize)) || (obj.lastParpoolSize ~= desiredNumWorkers)
                fprintf('Starting new parpool with %d workers\n', desiredNumWorkers);
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