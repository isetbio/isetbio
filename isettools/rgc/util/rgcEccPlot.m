% Plot RF diameter as a function of TEE
% James Golden
% 9/1/2015

clear

% pecc = [2.5 7.5 15 25 35];
% pcrad = [.03 .05 .07 .09 .15];
% pcrad_um = 225.3*pcrad;
% psrad = [.18 .43 .54 .73 .65];
% 
% figure; scatter(pecc, pcrad_um);
% 
% mecc = [5 15 25];
% mcrad = [.10 .18 .23];
% msrad = [.72 1.19 .58];

ecc_ch = [0.5 10];
dia_ch = [25  275];
m = (dia_ch(2)-dia_ch(1))/(ecc_ch(2)-ecc_ch(1));
yint = 25 - m*0.5;

xecc = 0:0.2:15; 
diaest = m*xecc + yint;

theta = 0:0.1:2*pi+0.1;
rad = 0:0.2:15;

theta_gr = repmat(theta,1,length(rad));
rad_gr = repmat(rad,1,length(theta));

figure; polar(theta_gr,rad_gr,'x');

rad_gr_eff = rad_gr;

ellip_pts = find(theta_gr>(pi/2) & theta_gr<(3*pi/2));

[xrad, yrad] = pol2cart(theta_gr, rad_gr);

aspect_ratio = 0.61
rad_gr_eff(ellip_pts) = sqrt((xrad(ellip_pts)*aspect_ratio).^2 + yrad(ellip_pts).^2);

for i = 1:length(theta_gr)
    tee(i) = retinalLocationToTEE((180/pi)*theta_gr(i), rad_gr(i), 1);
end

% figure; polar(theta_gr,rad_gr_eff,'x');

xradrs = reshape(xrad, length(theta), length(rad));
yradrs = reshape(yrad, length(theta), length(rad));
rad_gr_eff_rs = reshape(rad_gr_eff, length(theta), length(rad));
figure; contourf(xradrs, yradrs,rad_gr_eff_rs, 16);


tee_rs = reshape(tee, length(theta), length(rad));
figure; contourf(xradrs, yradrs, tee_rs, 16);

