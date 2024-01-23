function cal = CalibrateAmbDrvr(cal,USERPROMPT,whichMeterType,blankOtherScreen)
% cal =  CalibrateAmbDrvr(cal,USERPROMPT,whichMeterType,blankOtherScreen)
%
% This script does the work for monitor ambient calibration.

% 4/4/94    dhb     Wrote it.
% 8/5/94    dhb,ccc More flexible interface.
% 9/4/94    dhb     Small changes.
% 10/20/94  dhb     Add bgColor variable.
% 12/9/94   ccc     Nine-bit modification
% 1/23/95   dhb     Pulled out working code to be called from elsewhere.
%           dhb     Make user prompting optional.
% 1/24/95   dhb     Get filename right.
% 12/17/96  dhb, jmk  Remove big bug.  Ambient wasn't getting set.
% 4/12/97   dhb     Update for new toolbox.
% 8/21/97   dhb     Don't save files here.
%                   Always measure.
% 4/7/99    dhb   NINEBIT -> NBITS
%           dhb   Handle noMeterAvail, RADIUS switches.
% 9/22/99   dhb, mdr  Make boxRect depend on boxSize, defined up one level.
% 12/2/99   dhb   Put background on after white box for aiming.
% 8/14/00   dhb   Call to CMETER('Frequency') only for OS9.
% 8/20/00   dhb   Remove bits arg to SetColor.
% 8/21/00   dhb   Remove RADIUS arg to MeasMonSpd.
% 9/11/00   dhb   Remove syncMode code, any direct refs to CMETER.
% 9/14/00   dhb   Use OpenWindow to open.
%           dhb   Made it a function.
% 7/9/02    dhb   Get rid of OpenWindow, CloseWindow.
% 9/23/02   dhb, jmh  Force background to zero when measurements come on.
% 2/26/03   dhb   Tidy comments.
% 4/1/03    dhb   Fix ambient averaging.
% 8/19/12   dhb   Add codelet suggested by David Jones to clean up at end.  See comment in CalibrateMonSpd.
% 8/19/12   mk    Rewrite setup and clut code to be able to better cope with all
%                 the broken operating systems / drivers / gpus and to also
%                 support DataPixx/ViewPixx devices.
% 6/30/23   mk    Use new clut mapping to fix this mess on standard gpus.

global g_usebitspp;

% If the global flag for using Bits++ is empty, then it hasn't been
% initialized and default it to 0.
if isempty(g_usebitspp)
    g_usebitspp = 0;
end

% User prompt
if USERPROMPT
    if cal.describe.whichScreen == 0
        fprintf('Hit any key to proceed past this message and display a box.\n');
        fprintf('Focus radiometer on the displayed box.\n');
        fprintf('Once meter is set up, hit any key - you will get %g seconds\n',...
                cal.describe.leaveRoomTime);
        fprintf('to leave room.\n');
        KbStrokeWait(-1);
    else
        fprintf('Focus radiometer on the displayed box.\n');
        fprintf('Once meter is set up, hit any key - you will get %g seconds\n',...
                cal.describe.leaveRoomTime);
        fprintf('to leave room.\n');
    end
end

% Blank other screen, if requested:
if blankOtherScreen
    % We simply open an onscreen window with black background color:
    Screen('OpenWindow', cal.describe.whichBlankScreen, 0);
end

% Setup screen to be measured
% ---------------------------

% Prepare imaging pipeline for Bits+ Bits++ CLUT mode, or DataPixx/ViewPixx
% L48 CLUT mode (which is pretty much the same). If such a special output
% device is used, the Screen('LoadNormalizedGammatable', win, clut, 2);
% command uploads 'clut's into the device at next Screen('Flip'), taking
% care of possible graphics driver bugs and other quirks:
PsychImaging('PrepareConfiguration');

if g_usebitspp == 0
    % Setup for imaging pipeline based clut mapping:
    PsychImaging('AddTask', 'AllViews', 'EnableCLUTMapping', 256, 1);
end

if g_usebitspp == 1
    % Setup for Bits++ CLUT mode. This will automatically load proper
    % identity gamma tables into the graphics hardware and into the Bits+:
    PsychImaging('AddTask', 'General', 'EnableBits++Bits++Output');
end

if g_usebitspp == 2
    % Setup for DataPixx/ViewPixx etc. L48 CLUT mode. This will
    % automatically load proper identity gamma tables into the graphics
    % hardware and into the device:
    PsychImaging('AddTask', 'General', 'EnableDataPixxL48Output');
end

% Open the window:
[window, screenRect] = PsychImaging('OpenWindow', cal.describe.whichScreen, 0);

% Only need to set identity mapping for standard display mode, as all needed setup
% is already done automatically for CRS, VPixx devices etc. by PsychImaging():
if g_usebitspp == 0
    % Upload an identity gamma ramp into all channels. Upload a LUT with exactly
    % as many slots as the operating system + gpu driver expects, so we pass through
    % the full framebuffer precision one-to-one, as much as framebuffer depths allows.
    % Don't use LoadIdentityClut here, as it configures for a 8 bpc identity lut,
    % suitable for special display devices, but may lose a bit of precision for
    % conventional display setups:
    [~, ~, nslots] = Screen('ReadNormalizedGammaTable', window);
    Screen('LoadNormalizedGammaTable', window, (linspace(0, 1, nslots)' * ones(1, 3)), 0);
end

if (cal.describe.whichScreen == 0)
    HideCursor(window);
end

% Load zero theClut into device:
theClut = zeros(256,3);
Screen('LoadNormalizedGammaTable', window, theClut, 2);
Screen('Flip', window);

% Draw a box in the center of the screen
boxRect = [0 0 cal.describe.boxSize cal.describe.boxSize];
boxRect = CenterRect(boxRect, screenRect);
theClut(2,:) = [1 1 1];
Screen('FillRect', window, 1, boxRect);
Screen('LoadNormalizedGammaTable', window, theClut, 2);
Screen('Flip', window, 0, 1);

% Wait for user
if USERPROMPT == 1
    KbStrokeWait(-1);
    fprintf('Pausing for %d seconds ...', cal.describe.leaveRoomTime);
    WaitSecs(cal.describe.leaveRoomTime);
    fprintf(' done\n');
end

% Put in appropriate background.
theClut(1,:) = cal.bgColor';
Screen('LoadNormalizedGammaTable', window, theClut, 2);
Screen('Flip', window, 0, 1);

% Start timing
t0 = clock;

ambient = zeros(cal.describe.S(3), 1);
for a = 1:cal.describe.nAverage
    % Measure ambient
    ambient = ambient + MeasMonSpd(window, [0 0 0]', cal.describe.S, 'on', whichMeterType, theClut);
end
ambient = ambient / cal.describe.nAverage;

% Close the screen, restore cluts:
if g_usebitspp
    % Load identity clut on Bits++ / DataPixx et al.:
    BitsPlusPlus('LoadIdentityClut', window);
    Screen('Flip', window);
end

% Restore graphics card gamma tables to original state, close windows, etc.
sca;

% Report time:
t1 = clock;
fprintf('CalibrateAmbDrvr measurements took %g minutes\n', etime(t1,t0)/60);

% Update structure
Smon = cal.describe.S;
Tmon = WlsToT(Smon);
cal.P_ambient = ambient;
cal.T_ambient = Tmon;
cal.S_ambient = Smon;

% Done:
return;
