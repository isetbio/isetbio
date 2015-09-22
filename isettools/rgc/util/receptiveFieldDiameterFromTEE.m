function receptiveFieldDiameterParasol2STD = receptiveFieldDiameterFromTEE(tee)
% receptiveFieldDiameter2STD converts the temporal equivalent eccentricity
% (TEE) to a value for the RF diameter in um at the specified TEE. The RF
% diameter can then be used to bulid RGC RFs for


% NEED TO DO ACTUAL REGRESSION
% VERY ROUGH ESTIMATE FROM CHICHILNISKY & KALMAR (2002) SEE FIG. 5 PAGE 2741
% ON & OFF PARASOL CELLS 

scaleFactor = 1.57; % DF diameter = 1.57*(DF diameter)

ecc = [0.5 10]; dia2STD = [25  275]/scaleFactor;
m = (dia2STD(2)-dia2STD(1))/(ecc(2)-ecc(1));
yint = dia2STD(1) - m*ecc(1);

receptiveFieldDiameterParasol2STD = m*tee + yint;