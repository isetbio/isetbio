function [fractionBleached] = ComputePhotopigmentBleaching(irradiance,receptortype,units,source,initialFraction,timeUnits)
% [fractionBleached] = ComputePhtopigmentBleaching(irradiance,[receptorType],[units],[source],[initialFraction],[timeUnits])
%
% Compute fraction of photopigment bleached, given irradiance of light
% reaching the eye.
%
% There are two distinct uses, controlled by the value of initialFraction.
%
% Usage 1 - If initialFraction is not passed or is empty, the steady state
% fraction of pigment bleached is returned for each irradiance in the
% passed input irradiance.
%
% Usage 2 - If initialFraction is passed as a scalar, this is taken as the
% time zero fraction bleached, and inputirradiance is taken to be the time
% variying irradiance, with fractionBleached the time varying fraction
% bleached.
%
% When time varying signals are handled, the unit of time is as specified
% by the timeunits argument (default, msec).
%
% As far as I can tell, the fundemantal measurements of the half-bleach
% constant for human cones were made by Rushton and Henry (1968, Vision Research,
% 8, 617-631). This fact I learned from CVRL
% (http://www.cvrl.org/database/text/intros/introbleaches.htm).
%
% I am pretty sure that the Rushton and Henry measurements were made for
% 560 nm light, and they give (see their Figure 2) a half-bleach constant
% of 4.3 log10 trolands (20,000 td). This number is also given in Boynton
% and Kaiser, Human Color Vision, 2nd edition, pp 211 and following.  This
% is the number you get here if you specify 'Boynton' as the source for the
% cone bleaching data.
%
% Elsewhere in the same paper, Rushton and Henry use another method and
% come up with 29167 td) as the half-bleach constant.  You can get this
% number by specifying 'RushtonHenryAlt' as the source for the cone
% constants.  If you use isomemerizations here, the constant is scaled up
% from the one derived above by the ratio (29167/(10^4.3)).
%
% It's probably fine to compute bleaching for L and M cones given retinal
% illuminance in trolands, given that these are effects that matter over
% log10 units.  But trolands are not going to help much for the S-cones.
% According to CVRL there aren't good measurements for the half-bleaching
% constant for S cones because putting enough short-wavelength light onto
% the retina to bleach the S cones is not good for the eyes.
%
% None-the-less, it seems nice to have this routine written so that it will
% return a number if you give it irradiance either in trolands or in
% isomerizations/cone-sec.  For 560 nm light and the CIE 10 deg
% fundamentals, I compute that 1 td is 137 isomerizations/cone-sec for L
% cones and 110 isomerizations/cone-sec for M cones.  Take the weighted
% average value of (2*L + 1*M)/3 = 128 and multiply by (10.^4.3) to get a
% half-bleach constant in isomerizations/cone-sec of  2.55e+06 (6.4 log10
% isomerizations/cone-sec). [Computations done 6/2/14 using
% IsomerizationsInEyeDemo and setting the fundamentals to 'CIE10deg' and
% wavelength to 560 nm by hand in the code.  These are for the 'Boynton'
% source.]
%
% [ASIDE: I used 10 deg fundamentals to compute the bleaching constant
% expressed in terms of isomerizations, because I figure that Rushton's
% measurements are based on a fairly large field.  Because the macular
% pigment absorbs a fair amount of light, this matters.  If I compute
% instead with 2-deg fundamentals, I get that 1 td is 23.7 L cone
% isomerizations/cone-sec and 19.5 M cone isomerizations/cone-sec.   These
% two numbers are ballpark consistent with Rodiek page 475 who gives 18.3
% and 15.9 for a monochromatic 540 THz light (555 nm)].
%
% This paper
%   Burkhardt, D. A. "Light adaptation and photopigment
%   bleaching in cone photoreceptors in situ in the retina
%   of the turtle." Journal of Neuroscience 14.3 (1994):
%   1091-1105.
% provides a half bleach constant for turtle cones of 5.57 expressed
% in log10 R*/um2/sec, which could with some work be converted to
% isomerizations/cone/sec for turtle cones. But it's not
% clear you want to use that number unless you are studying turtle.
%
% This routine will do the computation either on the basis of input in
% trolands or input in isomerization/cone-sec, using the appropirate
% constant as above.  Note that the computation of isomerizations takes
% into account lens and macular pigment, while the troland value is the
% straight troland value.  A second advantage of using units of
% isomerizations/cone-sec is that you can compute this for other regions of
% the visual field and presumably the numbers will be about right.  You can
% also compute for S-cones on the assumption that the half-bleach constant
% is the same for S-cones as for L- and M- cones.
%
% As far as I can tell, the computations and analysis of bleaching do not
% take into account changes in isomerization rate that occur because of
% change in spectral sensitivity of cones with bleaching.  That is, the
% measurements are simply of steady state pigment density and are modeled
% with a formula that assumes monochromatic light (see treatment in
% Boynton).
%
% irradiance    -- retinal irradiance specified as determined by units. If
%                  initialFraction is empty, this is a single number and
%                  steady state bleaching fraction is returned.  If
%                  initialFraction is a number, then this is a time series
%                  of irradiance versus time, and fraction bleached for the
%                  same times is returned.
%
% receptortype
%   'cones'     -- computations for cones. [Default]
%   'rods'      -- computations for rods.
%
% units         -- units of irradiance
%   'trolands'     input irradiance trolands.  Note that the computation
%                  only makes sense for L and M cones if this is the input.
%                  This is photopic trolands if receptor type is 'cones'.
%                  [Default].  It is scotopic trolands if receptor type is
%                  rods.
%   'isomerizations'  nominal isomerization rate in
%                  isomerizations/cone-sec, comptued taking into account
%                  pre-retinal absorption as well as nominal cone axial
%                  density.  But not taking into account any pigment
%                  bleaching.
%
% source        -- source of underlying data
%   'Boynton'      Boynton and Kaiser, Human Color Vision, 2nd edition,
%                  pp. 211 and following. See intro text above. [Default]
%   'RushtonHenryAlt' - Rushton and Henry's other half-bleach constant.
%                  Only available for 'cones'.
%   'WyszeckiStiles' - Wyszecki and Stiles (1982) give parameters on page
%
% initialFraction -- fraction of input bleached at time zero. If
%                 empty, steady state fraction bleached is
%                 returned. Default is empty.
%
% timeUnits     -- units for time step on input and output
%   'sec'          seconds
%   'tenth'        tenths of seconds
%   'hundredth'    hundredths of seconds
%   'msec'         millseconds [Default]
%   'tenthmsec'     tenths of milliseconds
%
% 05/23/14 dhb  Wrote it.
% 05/26/14 dhb  Clean up.
% 06/02/14 dhb  Take isomerizations number based on 2:1 L:M assumed ratio.
% 12/18/18 dhb  Modify header comments for possibility of passing time
%               varying signal.  This breaks old usage that allowed
%               computing steady state bleaching for a set of vector
%               inputs, but I think that is OK.
% 08/19/19 dhb  Added some information about Burkhardt (1994) to header
%               comment, and inserted a stub to use that information if
%               someone does the work to put the number from it into the
%               right units.
%          dhb  Reorganized some code relative to source switch statement. Because
%               there was only one case this didn't matter, but now I think it
%               is right if more cases.
% 01/09/21 dhb  Finish off adding kinetics.
% 01/12/21 dhb  Add in the alternate Rushton & Henry half-bleach constant
%               for cones.

% Examples:
%{
    % Running this example should produce a plot that
    % looks like Figure 6.3 (p. 212) in Boyton and Kaiser
    % (red curve in plot with blue overlay here).
    clear; close all;
    isomerizationsPerTroland = 128;
    trolands = logspace(0,7,1000);
    isomerizations = isomerizationsPerTroland*trolands;
    for ii = 1:length(trolands)
        fractionBleached(ii) = ComputePhotopigmentBleaching(trolands(ii),'cones','trolands','Boynton');
        fractionBleached1(ii) = ComputePhotopigmentBleaching(isomerizations(ii),'cones','isomerizations','Boynton');
    end
    figure; clf; hold on;
    plot(log10(trolands),fractionBleached,'r','LineWidth',6);
    plot(log10(trolands),fractionBleached1,'b','LineWidth',4);
%}
%{
    % Running this example should produce a pigment recovery plot that
    % with a red curves looks like Figure 6.2 in Boynton (first edition),
    % Note that the righthand scale for that figure is flipped with 0 
    % at the top and 1 at the bottom.  We plot fraction bleached rather
    % than fraction unbleached, so our plot looks like this but has scale
    % in the conventional order. The red curve is for rods.  The blue
    % is for cones and is faster as expected.
    clear; close all;
    timeSec = (1:60*40) - 1;
    timeMinutes = timeSec/(60);
    irradiance = zeros(60*40,1);
    initialFractionBleached = 1;
    fractionBleached = ComputePhotopigmentBleaching(irradiance,'rods','trolands','WyszeckiStiles',initialFractionBleached,'sec');
    fractionBleached1 = ComputePhotopigmentBleaching(irradiance,'cones','isomerizations','Boynton',initialFractionBleached,'sec');
    figure; clf; hold on;
    plot(timeMinutes,fractionBleached,'r','LineWidth',6);
    plot(timeMinutes,fractionBleached1,'b','LineWidth',6);
    xlabel('Time (min)');  
    ylabel('Fraction bleached');
%}
%{
    % This example shows that the steady state calculation (red) matches the
    % asymptote of the dynamic calculation (green and blue), and that the asymptote 
    % is independent of the starting fraction bleached.
    clear; close all;
    timeSec = (1:60*10) - 1;
    timeMinutes = timeSec/(60);
    trolands = 10^5;
    irradiance = trolands*ones(60*10,1);
    fractionBleachedSteady = ComputePhotopigmentBleaching(trolands,'cones','trolands','Boynton');
    initialFractionBleached = 0;
    fractionBleached0 = ComputePhotopigmentBleaching(irradiance,'cones','trolands','Boynton',initialFractionBleached,'sec');
    initialFractionBleached = 1;
    fractionBleached1 = ComputePhotopigmentBleaching(irradiance,'cones','trolands','Boynton',initialFractionBleached,'sec');
    figure; clf; hold on;
    plot(timeMinutes,fractionBleachedSteady*ones(size(timeMinutes)),'r','LineWidth',6);
    plot(timeMinutes,fractionBleached0,'g','LineWidth',4);
    plot(timeMinutes,fractionBleached1,'b','LineWidth',4);
    xlabel('Time (min)');  
    ylabel('Fraction bleached');
%}

%% Specify receptor type
if (nargin < 2 || isempty(receptortype))
    receptortype = 'cones';
end

%% Specify units
if (nargin < 3 || isempty(units))
    units = 'trolands';
end

%% Specify source
if (nargin < 4 || isempty(source))
    source = 'Boynton';
end

%% Specify initial fraction
if (nargin < 5 || isempty(initialFraction))
    initialFraction = [];
end

%% Time units, relevant if initialFraction is not empty
if (nargin < 6 || isempty(timeUnits))
    timeUnits = 'msec';
end

%% Do it
switch (receptortype)
    case 'cones'
        switch (source)
            case 'Boynton'
                switch (units)
                    case 'trolands'
                        Izero = 10^4.3;   % In Td, stimulus intensity that bleaches at the rate of 1/N
                        N = 120;          % Recovery time constant (seconds)
                    case 'isomerizations'
                        Izero = 10^6.4;   % In L/M cone isomerizations/sec
                        N = 120;          % Recovery time constant (seconds)
                    otherwise
                        error('Unkown input units specified');
                end

            case 'WyszeckiStiles'
                switch (units)
                    case 'trolands'
                        Q = 5e6;          % Photopic td sec
                        N = 130;          % Recovery time constant (seconds)
                        Izero = Q/N;      % Photopic td, stimulus intensity that bleaches at the rate of 1/N
                    otherwise
                        error('Unkown input units specified');
                end

            case 'RushtonHenryAlt'
                switch (units)
                    case 'trolands'
                        Izero = 29167;    % In Td, stimulus intensity that bleaches at the rate of 1/N
                        N = 120;          % Recovery time constant (seconds)
                    case 'isomerizations'
                        Izero = (29167/(10^4.3))*10^6.4; % In L/M cone isomerizations/sec
                        N = 120;          % Recovery time constant (seconds)
                    otherwise
                        error('Unkown input units specified');
                end

            case 'Burkhardt'
                switch (units)
                    case 'trolands'
                        error('Only know how to deal units in isomserization/sec for Burkhardt case');
                    case 'isomerizations'
                        error('You need to do some more work to put the Burkhardt number into the right units');
                        %Izero = 10^5.57;
                    otherwise
                        error('Unkown input units for cones specified');
                end

            otherwise
                error('Unknown source for cones specified');
        end

    case 'rods'
        switch (source)
            case 'Boynton'
                error('No Boynton units entered for rods yet');
                switch (units)
                    % Values of Izero below are set to 1 just as a
                    % placeholder.  Need to be filled in with correct
                    % values to do full set of rod computations.
                    case 'trolands'
                        Izero = 1;        % Should be in scotopic Td, stimulus intensity that bleaches at the rate of 1/N
                        N = 400;          % Recovery time constant (seconds): 400=rhodopsin, 120=L/M cones.
                    case 'isomerizations'
                        Izero = 1;        % Should be in isomerizations/sec
                        N = 400;          % Recovery time constant (seconds): 400=rhodopsin, 120=L/M cones.
                    otherwise
                        error('Unsupported input units for rods specified');
                end

            case 'WyszeckiStiles'
                switch (units)
                    % Values of Izero below are set to 1 just as a
                    % placeholder.  Need to be filled in with correct
                    % values to do full set of rod computations.
                    case 'trolands'
                        Q = 1.57e7;       % Scotopic Td-sec
                        N = 519;          % Recovery time constant (seconds)
                        Izero = Q/N;
                    otherwise
                        error('Unsupported input units for rods specified');
                end

            otherwise
                error('Unknown source for rods specified');
        end
    otherwise
        error('Unknown receptor type specified');
end

% Steady state calculation
if (isempty(initialFraction))
    if (length(irradiance) ~= 1)
        error('Irradiance should be passed as a scalar for steady state computation');
    end

    % Simple formula.  This formula appears to have a typo in the first
    % edition (Eqn. 6.5 is missing the division sign).
    fractionBleached = (irradiance./(irradiance + Izero));

    % Time varying calculation.
else
    % Take time resolution into account. Our underlying
    % computations are in seconds, at the specified time
    % resolution.
    switch (timeUnits)
        case 'sec'
            timeStep = 1;
        case 'tenth'
            timeStep = 0.1;
        case 'hundredth'
            timeStep = 0.01;
        case 'msec'
            timeStep = 0.001;
        case 'tenthmsec'
            timeStep = 0.0001;
        otherwise
            error('Unknown time units specified');
    end

    % Inline functions to make our lives easier
    %
    % This is Eqn. 6.5 from Boynton Human Color Vision
    % (first edition) and gives the change in fraction
    % unbleached (p) over time.
    dp_dt = @(p, I) (( (1-p)./N ) - ( (I.*p)./ (N.*Izero) ));

    % Set up fraction bleached over simulation time.  Note conversion
    % of passed fraction bleached to fraction unbleached form here.
    fractionUnbleached = zeros(size(irradiance));
    fractionUnbleached(1) = 1-initialFraction;

    for t=2:length(irradiance)
        fractionUnbleached(t) = fractionUnbleached(t-1) + ...
            dp_dt(fractionUnbleached(t-1),irradiance(t-1))*timeStep;
    end

    % Convert back to fraction bleached
    fractionBleached = 1 - fractionUnbleached;

end
end




