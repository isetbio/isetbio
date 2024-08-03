function same = DescribeMonCal(calOrCalStruct,file,whichScreen)
% same = DescribeMonCal(calOrCalStruct,[file],[whichScreen])
% 
% Print descriptive information about a calibration 
% to the command window or file.
%
% Argument file is a standard Matlab file descriptor,
% see fopen.  If file arg is omitted or empty, printout
% goes to command window.
%
% If argument whichScreen is passed, a description of
% the current hardware is also printed.  In this case,
% returned boolean same indicates whether the calibration
% is consistent with the current hardware.  Boolean
% same is empty if whichScreen is not provided.
%
% 8/25/97  dhb, pbe  Wrote it.
% 7/3/98   dhb, pbe  Updated for cal.describe.
% 12/3/99  dhb, mpr  Fix check for calibration desription field.
% 8/1800   dhb       Add whichScreen arg, same return.
% 6/29/02  dgp       Use new version of Screen VideoCard.
% 9/23/02  dhb, jms  Fix small bug in way driver is compared, presumably introduced 6/29/02.
% 9/29/08  dhb, tyl, ijk Update for OS/X, current computer stuff.
%                    Comparison of computer name skipped, because it seems to vary with login. 
% 6/24/11  dhb       Dump out gamma fit type and exponents if gamma function was fit with a simple power function.
% 5/28/13  dhb       Change output printed format to make it easier to paste into Doku wiki.
% 5/08/14  npc       Modifications for accessing calibration data using a @CalStruct object.
%                    The first input argument can be either a @CalStruct object (new style), or a cal structure (old style).
%                    Passing a @CalStruct object is the preferred way because it results in 
%                    (a) less overhead (@CalStruct objects are passed by reference, not by value), and
%                    (b) better control over how the calibration data are accessed.

% Default args
if (nargin < 2 || isempty(file))
	file = 1;
end
if (nargin < 3 || isempty(whichScreen))
	file = 1;
	whichScreen = [];
end
same = [];

% Specify @CalStruct object that will handle all access to the calibration data.
[calStructOBJ, inputArgIsACalStructOBJ] = ObjectToHandleCalOrCalStruct(calOrCalStruct);
if (~inputArgIsACalStructOBJ)
     % The input (calOrCalStruct) is a cal struct. Clear it to avoid  confusion.
    clear 'calOrCalStruct';
end
% From this point onward, all access to the calibration data is accomplised via the calStructOBJ.

fprintf('Calibration:\n');
fprintf(file,'  * Computer: %s\n',                      calStructOBJ.get('computer'));         % cal.describe.computer;
fprintf(file,'  * Screen: %d\n',                        calStructOBJ.get('whichScreen'));      %cal.describe.whichScreen);
fprintf(file,'  * Monitor: %s\n',                       calStructOBJ.get('monitor'));          %cal.describe.monitor);
fprintf(file,'  * Video driver: %s\n',                  calStructOBJ.get('driver'));           %cal.describe.driver);
fprintf(file,'  * Dac size: %g\n',                      calStructOBJ.get('dacsize'));          %cal.describe.dacsize);
fprintf(file,'  * Frame rate: %g hz\n',                 calStructOBJ.get('hz'));               %cal.describe.hz);
fprintf(file,'  * Calibration performed by %s\n',       calStructOBJ.get('who'));              %cal.describe.who);
fprintf(file,'  * Calibration performed on %s\n',       calStructOBJ.get('date'));             %cal.describe.date);
fprintf(file,'  * Calibration program: %s\n',           calStructOBJ.get('program'));          %cal.describe.program);
fprintf(file,'  * Comment: %s\n',                       calStructOBJ.get('comment'));          %cal.describe.comment);
fprintf(file,'  * Calibrated device has %g primaries\n',calStructOBJ.get('nDevices'));         %cal.nDevices);
fprintf(file,'  * Gamma fit type %s\n',                 calStructOBJ.get('gamma.fitType'));    %cal.describe.gamma.fitType);


if (strcmp(calStructOBJ.get('gamma.fitType'),'simplePower'))
    exponents = calStructOBJ.get('gamma.exponents');
    fprintf(file,'  * Simple power gamma exponents are: %0.2f, %0.2f, %0.2f\n',...
        exponents(1), exponents(2), exponents(3));
end
fprintf(file,'\n');

% Current configuration
if (~isempty(whichScreen))     
    computerInfo = Screen('Computer');
    computer     = sprintf('%s''s %s, %s', computerInfo.consoleUserName, computerInfo.machineName, computerInfo.system);
    driver       = sprintf('%s %s','unknown_driver','unknown_driver_version');
    dacsize      = ScreenDacBits(whichScreen);
    hz           = Screen('NominalFrameRate',whichScreen);
    
    same = 1;
    fprintf('Current configuration:\n');
    fprintf(file,'  * Computer: %s\n',computer);
    if (~streq(computer,calStructOBJ.get('computer')))
        %same = 0;
    end
    fprintf(file,'  * Screen: %d\n',whichScreen);
    if (whichScreen ~= calStructOBJ.get('whichScreen'))
        save = 0;
    end
    fprintf(file,'  * Video driver: %s\n',driver);
    fprintf(file,'  * Dac size: %g\n',dacsize);
    if (dacsize ~= calStructOBJ.get('dacsize'))
        same = 0;
    end
    fprintf(file,'  * Frame rate: %g hz\n',hz);
    if (abs(hz - calStructOBJ.get('hz')) > 0.5)
        same = 0;
    end
end
