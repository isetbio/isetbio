function navarroAccomm = convertToNavarroAccomm(inputAccomm)
% Convert accommodations to use Navarro with data from other models
%
% Syntax:
%   navarroAccomm = convertToNavarroAccomm(inputAccomm)
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
%    inputAccomm   - Numeric. The non-Navarro accommodation. This should be
%                    within a reasonable range for accommodation, aka it
%                    should be between 0.44 and 9.25.
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

% Catch the special case of 0 accommodation. 
% if(inputAccomm == 0)
%     navarroAccomm = 0;
%     return;
% end

% We can't interpolate past ~0.4 diopters and over ~9.25. 
if(inputAccomm < 0.44) 
    %navarroAccomm = inputAccomm;
    navarroAccomm = 0;
    return;
elseif(inputAccomm > 9.25)
    warning(['At the moment, we don''t have data for accommodation ', ...
    'greater than 9.25 dpt. This may change in the near future, but ', ...
    'for now we set accommodation to 9.25 dpt.']);
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
