function isSame = CompareMonCal(calOrCalStruct1,calOrCalStruct2,IGNOREDATE)
% isSame = CompareMonCal(calOrCalStruct1,calOrCalStruct2,[IGNOREDATE])
%
% Checks if the two calibrations are the same.  Useful
% for preventing blunders if you have programs that
% precompute and save quantities based on monitor calibrations.
% In that case, this can be used to ensure that current 
% calibration matches the one used to do the pre-computing.
%
% Checks date/time, screen, and computer.  Could check the
% actual data, but that seems like overkill.
% 
% 9/17/97  pbe       Wrote it. 
% 9/18/97  pbe, dhb  Modify interface, change name.
% 1/16/98  dhb       Add any around string compares, necessary for desired effect.
% 1/21/98  dhb       Add IGNOREDATE flag.
% 3/10/98  dhb	     Change name to CompareMonCal.
% 7/3/98   dhb, pbe  Change for cal.describe format.
% 5/28/14  npc       Modifications for accessing calibration data using a @CalStruct object.
%                    The first two input arguments can be either a @CalStruct object (new style), or a cal structure (old style).
%                    Passing @CalStruct objects is the preferred way because it results in 
%                    (a) less overhead (@CalStruct objects are passed by reference, not by value), and
%                    (b) better control over how the calibration data are accessed.

if (nargin < 3 || isempty(IGNOREDATE))
	IGNOREDATE = 0;
end

% Specify @CalStruct objects that will handle all access to the calibration data.
[calStruct1OBJ, inputArg1IsACalStructOBJ] = ObjectToHandleCalOrCalStruct(calOrCalStruct1);
if (~inputArg1IsACalStructOBJ)
     % The input (calOrCalStruct) is a cal struct. Clear it to avoid  confusion.
    clear 'calOrCalStruct1';
end
[calStruct2OBJ, inputArg2IsACalStructOBJ] = ObjectToHandleCalOrCalStruct(calOrCalStruct2);
if (~inputArg2IsACalStructOBJ)
     % The input (calOrCalStruct) is a cal struct. Clear it to avoid  confusion.
    clear 'calOrCalStruct2';
end
% From this point onward, all access to the calibration data is accomplised via the calStructOBJ.

% Get necessary calibration data
date1 = calStruct1OBJ.get('date');
date2 = calStruct2OBJ.get('date');

whichScreen1 = calStruct1OBJ.get('whichScreen');
whichScreen2 = calStruct2OBJ.get('whichScreen');

computer1 = calStruct1OBJ.get('computer');
computer2 = calStruct2OBJ.get('computer');

driver1 = calStruct1OBJ.get('driver');
driver2 = calStruct2OBJ.get('driver');

dacsize1 = calStruct1OBJ.get('dacsize');
dacsize2 = calStruct2OBJ.get('dacsize');

% Compare
isSame = 1;
if (~IGNOREDATE)
	if (~streq(date1, date2))
		%fprintf(1,'CompareCal:\n');
		%fprintf(1,'\tcal1 calibration date: %s',date1);
		%fprintf(1,'\tcal2 calibration date: %s',date2);
		isSame = 0;
	end
end
if (whichScreen1 ~= whichScreen2)
	%fprintf(1,'CompareCal:\n');
	%fprintf(1,'\tcal1 calibration screen: %g\n',whichScreen1);
	%fprintf(1,'\tcal2 calibration screen: %g\n',whichScreen2);
	isSame = 0;
end	
if (~streq(computer1,computer2))
	%fprintf(1,'CompareCal:\n');
	%fprintf(1,'\tcal1 computer: %s',computer1);
	%fprintf(1,'\tcal2 computer: %s',computer2);
	isSame = 0;
end
if (~streq(driver1,driver2))
	%fprintf(1,'CompareCal:\n');
	%fprintf(1,'\tcal1 driver: %s\n',driver1);
	%fprintf(1,'\tcal2 driver: %s\n',driver2);
	isSame = 0;
end
if (dacsize1 ~= dacsize2)
	%fprintf(1,'CompareCal:\n');
	%fprintf(1,'\tcal1 DAC size: %g\n',dacsize1);
	%fprintf(1,'\tcal2 DAC size: %g\n',dacsize2);
	isSame = 0;
end

