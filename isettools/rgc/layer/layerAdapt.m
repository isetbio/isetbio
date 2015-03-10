function adaptTS = layerAdapt(layer)
% Adaptation (and nonlinearity) function for the RGC voltage timeseries
%
%   adaptTS = layerAdapt(layer)
%
% We have to think through whether we want positive and negative voltages
% preserved at this point.  In the first implementation we had negative
% voltages but scaled the range.
%
% This is like assuming that the background is at 0 and we are coding
% differential voltages from the background.  That could be OK.
% 
% This function scales the time series so that 
%
%  * The max is equal to the voltage swing
%  * There is a slight nonlinearity (logistic)
%
% Over time, we will have to implement this as a dynamic process.
%
%(c) Stanford VISTA Team 2011

if notDefined('layer'), error('layer structure required'); end

adaptTS = layer.get('lin ts');
vSwing = layer.get('vSwing');

% This just rescales to vSwing. We should check the literature, think about
% saturation, all kinds of stuf
%  mx = max(adaptTS(:));
%  mn = min(adaptTS(:));
%  adaptTS = vSwing * (adaptTS - mn)/(mx - mn); 

% For a while, we just did this.
% Negative is still down and positive is still up
% But the range is shifted and scaled to 0,1.
adaptTS = ieScale(adaptTS,vSwing);
% max(adaptTS(:)), min(adaptTS(:))

% Then, we put this Michaelis-Menten nonlinearity in. 
% Should we add parameters for this function?
% The semi-saturation and an exponent could be parameters.

% Maybe this should be a pointer to a function.
% We are trying to make a saturating nonlinearity.  This could work.
%
% adaptTS = adaptTS - min(adaptTS(:));
% adaptTS = vSwing * (adaptTS ./ (adaptTS + 0.5*vSwing));
% adaptTS = adaptTS + min(adaptTS(:));
%
% v = (0.0005:0.001:vSwing);
% vcNewGraphWin; semilogx(v,vSwing * (v ./ (v + 0.05*vSwing)))

% ANother possible parameterization using a logistic
% v = (-vSwing:0.001:vSwing);
% vcNewGraphWin; semilogx(v,vSwing ./ (1 + exp(-(v*100))));

% adaptTS = vSwing ./ (1 + 2*exp(-adaptTS ));


return
