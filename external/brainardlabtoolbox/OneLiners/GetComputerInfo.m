function computerInfo = GetComputerInfo
% GetComputerInfo - Gets basic information about the computer and user.
%
% Syntax:
% computerInfo = GetComputerInfo
%
% Description:
% Gets the following computer/user information
% 1. CPU architecture, e.g. i386
% 2. Matlab platform, e.g. MACI64
% 3. Physical memory in MB.
% 4. User short name
% 5. User long name
% 6. OS version, e.g. 10.6.8
% 7. Network name, e.g. computer1.example.com
% 8. Localhost name, e.g. Computer1
% 9. MAC address of the 1st ethernet port.
%
% Output:
% computerInfo (struct) - A struct containing the computer/user data.

% Gets the type of Mac architecture we're on as reported by the computer.
[~, architecture] = system('uname -p');
computerInfo.architecture = architecture(1:end-1);

% Gets the Matlab platform: MAC (PPC), MACI (32-bit Intel), or MACI64
% (64-bit Intel).
computerInfo.MatlabPlatform = computer;

% Get the number of CPUs.
[~, numCPUs] = system('sysctl -n hw.ncpu');
computerInfo.numCPUs = str2double(numCPUs(1:end-1));

% Get the amount of physical memory in MB.
[~, physicalMemory] = system('sysctl -n hw.physmem');
computerInfo.physicalMemory = str2double(physicalMemory(1:end-1)) / 1024^2;

% Gets the console username (short name) of the current user.
[~, userShortName] = system('id -un');
computerInfo.userShortName = userShortName(1:end-1);

% Gets the long name of the user.
directoryCommand = sprintf('dscl . read /users/%s RealName', ...
	computerInfo.userShortName);
[~, userLongName] = system(directoryCommand);
computerInfo.userLongName = userLongName(12:end-1);

% Get the OS version.
[~, OSVersion] = system('sw_vers -productVersion');
computerInfo.OSVersion = OSVersion(1:end-1);

% Get the network name of the computer.
[~, networkName] = system('uname -n');
computerInfo.networkName = networkName(1:end-1);

% Get the local host name.
[~, localHostName] = system('scutil --get LocalHostName');
computerInfo.localHostName = localHostName(1:end-1);

% Get the MAC address.
[~, s] = system('ifconfig en0 ether');
i = strfind(s, 'ether ');
computerInfo.MACAddress = upper(s(i+length('ether '):end-2));
