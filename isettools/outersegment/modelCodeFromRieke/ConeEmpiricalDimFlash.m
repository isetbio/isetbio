function [fit] = ConeEmpiricalDimFlash(coef, t)
% function [fit] = ConeEmpiricalDimFlash(coeffs, time)
% Following Schnapf, Baylor, 1990: Damped oscillator with S-shaped onset
% ScFact = coef(1); %Scaling Factor
% TauR = coef(2); %Rising Phase Time Constant
% TauD = coef(3); %Dampening Time Constant
% TauP=coef(4); %Oscillator Period
% Phi=coef(5); %Oscillator Phase
% fit = ScFact .* (((t./TauR).^3)./(1+((t./TauR).^3))) .* exp(-((t./TauD))).*cos(((2.*pi.*t)./TauP)+(2*pi*Phi/360));
% Created Apr_2011 Angueyra
% Modified Apr_2011 Fred (replaced gaussian by exponential decay)


ScFact = coef(1);% Scaling Factor
TauR = coef(2); %Rising Phase Time Constant
TauD = coef(3); %Dampening Time Constant
TauP = coef(4); %Period
Phi = coef(5); %Phase


% fit = ScFact .* (((t./TauR).^3)./(1+((t./TauR).^3))) .*exp(-((t./TauD).^2)).*cos(((2.*pi.*t)./TauP)+(2*pi*Phi/360));

fit = ScFact .* (((t./TauR).^3)./(1+((t./TauR).^3))) .* exp(-((t./TauD))) .* cos(((2.*pi.*t)./TauP)+(2*pi*Phi/360));

% fit = ScFact .* (((t./TauR).^4)./(1+((t./TauR).^4))) .* exp(-((t./TauD))) .* cos(((2.*pi.*t)./TauP)+(2*pi*Phi/360));
