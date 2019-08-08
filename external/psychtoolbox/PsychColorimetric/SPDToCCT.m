function CCT = SPDToCCT(SPD,S_SPD)
% CCT = SPDToCCT(SPD,S_SPD)
%
% Calculate Correlated Colour Temperature (CCT) from a spectral power
% distribution.
%
% Works by computing 1960 uv coordinates of black body radiators and then
% finding closest to the uv coordinates of the passed spectral power
% distribution.
%
% The source code contains examples.  Edit and then select and execute.

% History:
%
%                Written by Danny Garside 
% 05/04/19  dhb  Try to fix line breaks, tune up comments.
% 05/10/19  dhb  Use bsxfun rather than -, because dimension expansion does
%                not work in older versions of Matlab.

% Examples:
%{ 
    % Demonstrate proper inversion on black body radiators
    S_SPD = [380,5,81];
    theCCT = 6500;
    SPD = GenerateBlackBody(theCCT,SToWls(S_SPD));
    foundCCT = SPDToCCT(SPD,S_SPD);
    if (theCCT == foundCCT)
        fprintf('Self-inversion works\n');
    else
        fprintf('Input CCT is %d, found is %d\n',theCCT,foundCCT);
    end
%}
%{
    % Demonstrate approximate inversion on CIE daylights
    % This comes out at 6501 rather than 6500 for unknown reasons.
    load B_cieday
    theCCT = 6500;
    SPD = GenerateCIEDay(theCCT,B_cieday);
    foundCCT = SPDToCCT(SPD,S_cieday);
    if (theCCT == foundCCT)
        fprintf('Self-inversion works\n');
    else
        fprintf('Input CCT is %d, found is %d\n',theCCT,foundCCT);
    end
%}
%{
    % Calculate a bunch
    load spd_houser.mat
    CCT = SPDToCCT(spd_houser,S_houser);
%}


%% Calculate colorimetry
load T_xyz1931.mat T_xyz1931 S_xyz1931
if any(S_xyz1931 ~= S_SPD)
    T_xyz1931 = SplineCmf(S_xyz1931,T_xyz1931,S_SPD); %set to same wavelength range and sampling interval
end
SPD_XYZ = T_xyz1931*SPD;
SPD_uv = XYZTouv(SPD_XYZ,'Compute1960');

%% Calculate look ups
range = 1000:10000;
lookup_uv = zeros(2,length(range));
for i=1:length(range)
    lookup_uv(:,i) = XYZTouv(T_xyz1931*GenerateBlackBody(range(i),SToWls(S_SPD)),'Compute1960');
end

%% Find closest
for i=1:size(SPD,2)
    [~,minloc(i)] = min(sqrt(sum((bsxfun(@minus,lookup_uv,SPD_uv(:,i))).^2)));
    CCT(i) = range(minloc(i));
end

%% Check
if any(CCT == range(1) | CCT == range(end))
    warning(['You have exceeded the boundary of the temps searched for CCT. Current range is ',num2str(range(1)),'K:',num2str(range(end)),'K'])
end
