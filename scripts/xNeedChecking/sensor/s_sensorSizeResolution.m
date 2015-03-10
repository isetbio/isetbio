% s_SensorSizeResolution
%
% For a particular size sensor (quarterinch, halfinch), we can calculate
% how many pixels fit.
%
% Copyright ImagEval Consultants, LLC, 2005.

% halfinch, quarterinch
sSize = sensorFormats('halfinch');    % Returns sensor size in meters
pSize = (0.8:0.2:3)*1e-6;             % Pixel size in meters

r = zeros(size(pSize));
c = zeros(size(pSize));
m = zeros(size(pSize));

for ii=1:length(pSize)
    r(ii) =  sSize(1) / pSize(ii);
    c(ii) =  sSize(2) / pSize(ii);
    m(ii) = ieN2MegaPixel(r(ii)*c(ii));
end

vcNewGraphWin
plot(pSize*1e6,m,'-o')
grid on
xlabel('Pixel size (um)')
ylabel('Megapixel')

%% End
