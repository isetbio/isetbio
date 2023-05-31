function navarroAccomm = convertToNavarroAccomm(inputAccomm)
% Convert accommodations to use Navarro with data from other models
%
% Syntax:
%   navarroAccomm = convertToNavarroAccomm(inputAccomm)
%
% Description (TL):
%    We need this function because the input accommodation of the
%    Navarro model does not precisely match the distance to the point
%    in the scene we wish to accommodate to for 550 nm.
%
%    For example, if we use the 0.4 diopter Navarro model, the point
%    of best focus for 550 nm would be at 1013 mm, according to Zemax
%    (which we trust). The user probably wanted the point of best
%    focus to be 1/0.4 diopter = 1500 mm. This discrepancy probably
%    happens because Navarro fitted the accommodation equations to
%    match 5 and 10 diopters, but not other accommodation states.
%    Either that, or Navarro is somehow defining accommodation
%    differently than we are here.
%
%    We (TL) used Zemax and a few sample points to find the Navarro
%    accommodation state that most closely matches the user's desired
%    point of best focus. BW isn't sure where that calculation is,
%    however.
%
% Inputs:
%    inputAccomm   - Numeric. The non-Navarro accommodation. We can
%                    accept a reasonable range for accommodation, from
%                    0.44 D (2.27 m) to 9.25 (10 cm).  Anything
%                    smaller than 0.44D is set to 0 (Inf) and anything
%                    closer than 10 cm is set to 10 cm.  Warnings are
%                    issued by this code when the values are changed.
%
% Outputs:
%    navarroAccomm - Numeric. The Navarro model accommodation.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: XXX - The values listed below were sampled manually in Zemax.
%      Eventually we will have macros in Zemax that can sample this
%      relationship more finely.]
%

% Examples:
%{
% This describes the conversion that TL implemented in the form of
% some graphs.
%

% Input distance
inD = 0.2:0.05:2;  % Meters
inA  = 1./inD;     % Input accommodation

% The conversion
outA = zeros(size(inA));
for ii=1:numel(inA)
 outA(ii) = convertToNavarroAccomm(inA(ii));
end

% Accommodation conversion TL implemented
ieNewGraphWin;
plot(inA,outA,'-ko');
grid on; identityLine;
xlabel('input acc (D)'); ylabel('output acc (D)');

% Distance conversion
ieNewGraphWin;
plot(inD, 1 ./ outA,'-ko');
grid on; identityLine;
xlabel('input dist (m)'); ylabel('output dist (m)');

%}

% Catch the special case of 0 accommodation. 
if(inputAccomm == 0)
    navarroAccomm = 0;
    return;
end

% These are limits on our interpolation
% Past ~0.4 diopters (2.27 meters) and over ~9.25 (10 cm)
if(inputAccomm < 0.44) 
    warning('Accommodation < 0.44 (> 2.27 m) is set to 0 (~Inf m)');
    % treated as infinitely far away.
    navarroAccomm = 0;
    return;
elseif(inputAccomm > 9.25)
    warning('Accommodation > 9.25D (~10 cm) set to 9.25D.');
    navarroAccomm = 9.25;
    return;
end

accom = [0; 0.1500; 0.2400;
    0.40; 0.6300; 1.0000;
    1.5; 2.0000; 2.5; 3.0000;
    3.5; 4; 4.5; 5.0000; 6; 6.5;
    7.0000; 8; 8.5;
    9; 9.6; 10.0000];

dist550 = [2275.4; 1518.4; 1280.950;
    1013.579; 795.72; 603.02;
    466.509; 385.92; 331.787; 292.238;
    263.988; 244.850; 217.098; 200.036; 172.586; 161.308
    148.158; 131.835; 125.213;
    119.859; 112.590; 108.021];

% Plot
%{
figure;
plot(dist550,accom,'ro-'); hold on;
dist = 50:2300;
plot(dist,1./(dist*10^-3),'-');
xlabel('Focal Distance')
ylabel('Accommodation');
legend('Navarro model','Ideal');
grid on;
set(findall(gcf,'-property','FontSize'),'FontSize',18)
set(findall(gcf,'-property','LineWidth'),'LineWidth',2)

figure;
plot(1./(dist550*10^-3),accom,'bo-'); hold on;
xlabel('Desired Accommodation')
ylabel('Navarro Accommodation');
grid on;
set(findall(gcf,'-property','FontSize'),'FontSize',18)
set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
axis([0 10 0 10])
% axis image;
%}

% Linearly interpolate
navarroAccomm = interp1(1 ./ (dist550 * 10 ^ -3), accom, inputAccomm);

end
