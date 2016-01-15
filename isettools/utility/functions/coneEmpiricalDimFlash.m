function response = coneEmpiricalDimFlash(coef, t)
% response = coneEmpiricalDimFlash(coeffs, time)
%
% Following Schnapf, Baylor, 1990: Damped oscillator with S-shaped onset
%
% ScFact = coef(1); %Scaling Factor
% TauR = coef(2); %Rising Phase Time Constant
% TauD = coef(3); %Dampening Time Constant
% TauP=coef(4); %Oscillator Period
% Phi=coef(5); %Oscillator Phase
% fit = ScFact .* (((t./TauR).^3)./(1+((t./TauR).^3))) .* exp(-((t./TauD))).*cos(((2.*pi.*t)./TauP)+(2*pi*Phi/360));
%
% Apr_2011  Angueyra   Created 
% Apr_2011  Rieke      Replaced gaussian by exponential decay
% 1/8/16    dhb        Rename, clean

%% GIve names to the coefficient vector entries
ScFact = coef(1);% Scaling Factor
TauR = coef(2); %Rising Phase Time Constant
TauD = coef(3); %Dampening Time Constant
TauP = coef(4); %Period
Phi = coef(5); %Phase

%% Compute the response
response = ScFact .* (((t./TauR).^3)./(1+((t./TauR).^3))) .* exp(-((t./TauD))) .* cos(((2.*pi.*t)./TauP)+(2*pi*Phi/360));

%% Some earlier versions now not used
% response = ScFact .* (((t./TauR).^3)./(1+((t./TauR).^3))) .* exp(-((t./TauD).^2)).*cos(((2.*pi.*t)./TauP)+(2*pi*Phi/360));
% response = ScFact .* (((t./TauR).^4)./(1+((t./TauR).^4))) .* exp(-((t./TauD))).*cos(((2.*pi.*t)./TauP)+(2*pi*Phi/360));
