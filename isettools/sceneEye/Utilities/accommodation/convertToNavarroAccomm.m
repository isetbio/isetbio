function [ navarroAccomm ] = convertToNavarroAccomm(inputAccomm)
% Convert accommodations to use Navarro with data from other models
%
% Syntax:
%   [navarroAccomm] = convertToNavarroAccomm(inputAccomm)
%
% Description:
%    The reason we need this function is because the input accommodation to
%    the Navarro model does not precisely match the distance to the point
%    in the scene we wish to accommodate to (for 550 nm).
%
%    For example, if we use the 0.4 diopter Navarro model, the point of
%    best focus for 550 nm would be at 1013 mm, according to Zemax.
%    However, the user probably wanted the point of best focus to be 1/0.4
%    diopter = 1500 mm. This discrepancy probably happens because Navarro
%    fitted the accommodation equations to match 5 and 10 diopters, but not
%    other accommodation states. Either that, or Navarro is somehow
%    defining accommodation differently than we are here.
%
%    Nevertheless, here we use a few sample points to find the Navarro
%    accommodation state that most closely matches the user's desired point
%    of best focus.
%
% Inputs:
%    inputAccomm   - The non-Navarro accommodation
%
% Outputs:
%    navarroAccomm - The Navarro model accommodation.
%
% Notes:
%    * [Note: XXX - The values listed below were sampled manually in Zemax.
%      Eventually we will have macros in Zemax that can sample this
%      relationship more finely.]
%

% Catch the special case of 0 accommodation. 
if(inputAccomm == 0)
    navarroAccomm = 0;
    return;
end

% We can't interpolate past ~0.4 diopters. At that point the difference is
% minimal, so let's just use the original value
if(inputAccomm < 0.5)
    navarroAccomm = inputAccomm;
    return;
end

accom = [0; 0.1500; 0.2400; 0.40; 0.6300; 1.0000; 1.5; 2.0000; 3.0000;
    5.0000; 7.0000; 10.0000];

dist550 = [2275.4; 1518.4; 1280.950; 1013.579; 795.72; 603.02; 466.509;
       385.92; 292.23; 198.72; 151.26; 98.886];

% Linearly interpolate
navarroAccomm = interp1(1 ./ (dist550 * 10 ^ -3), accom, inputAccomm);

end
