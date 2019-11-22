function gauss = gauss2(halfWidthY, supportY, halfWidthX, supportX)
% Construct a 2D Gaussian matrix. 
%
% Syntax:
%   gMatrix = gauss2(halfWidthY, supportY, halfWidthX, supportX)
%
% Description:
%    You should always call this function with an odd number for the
%    support so that the Gaussian is symmetric (centered). In that case the
%    calculation here ends up being the same as the calculation used by
%    Matlab in fspecial. There is a slight difference for an even support.
%
%    The bivariate Gaussian formula is 
% 
%       g2 = (1/(2*pi*sx*sy)) * exp(-(1/2)*(x/sx)^2 + (y/sy)^2)
%
% Inputs:
%    halfWidthY - Numeric. The half width half max point across the Y axis.
%    supportY   - Numeric. The spatial support across the Y axis.
%    halfWidthX - (Optional) Numeric. If not provided, assumed same as
%                 halfWidthY. The half width half max point across X axis.
%    supportX   - (Optional) Numeric. If not provided, assumed same as
%                 supportY. The spatial support across the X axis.
%
% Outputs
%    gauss      - Matrix. A matrix of supportY by supportX, representing
%                 the created Gaussian matrix.
%
% Optional key/value pairs:
%    None.
%
% References:
%    * http://en.wikipedia.org/wiki/Multivariate_normal_distribution#...
%      Bivariate_case
%

% History:
%    XX/XX/03       Copyright ImagEval Consultants, LLC, 2003.
%    05/28/19  JNM  Documentation pass

% Examples:
%{
    % These casese yield the same result when the support is odd.
    % They differ when the support is even
    gG = gauss2(5, 25, 5, 25);
    mesh(1:25, 1:25, gG)
    gF = fspecial('gauss', 25 , ieHwhm2SD(5));
    mesh(1:25, 1:25, gF)
    mesh(gG - gF)

    gG = gauss2(2, 25, 2, 25);
    mesh(1:25, 1:25, gG)
    gF = fspecial('gauss', 25, ieHwhm2SD(2));
    mesh(1:25, 1:25, gF)
    mesh(gG - gF)

    % [Note: JNM - 3/11 failed, change everything to 11 to match and does
    % not error on Data dimensions (error using mesh (line 71), Data
    % dimensions must agree).]
    %
    % When the support is even, the difference is small. It occurs due to
    % a slight shift in the curve that we introduce to center the support.
    gG = gauss2(3, 11, 3, 11);
    mesh(1:11, 1:11, gG)
    gF = fspecial('gauss', 11 , ieHwhm2SD(3));
    mesh(1:11, 1:11, gF)
    mesh(gG - gF)

    gG = gauss2(8, 50, 8, 50);
    mesh(1:50, 1:50, gG)
    gF = fspecial('gauss', 50 , ieHwhm2SD(8));
    mesh(1:50, 1:50, gF)
    mesh(gG - gF)

    % We differ slightly fspecial, sigh. But not much. 
%}

if notDefined('halfWidthX'), halfWidthX = halfWidthY; end
if notDefined('supportX'), supportX = supportY;   end

% Calculate the spatial support. This is the place where we have a slight
% difference with Mathworks. Probably, we should just do it their way.
x = (1:supportX) - round(supportX / 2);
y = (1:supportY) - round(supportY / 2);
[X, Y] = meshgrid(x, y);

% Calculate the standard deviations from the half widths. 
sdX = ieHwhm2SD(halfWidthX);
sdY = ieHwhm2SD(halfWidthY);

% Gaussian formula, finally, apart from the scale factor
gauss = exp(-(1 / 2) * ((X / sdX) .^ 2 + (Y / sdY) .^ 2));

% Normalize to unit area
gauss = gauss / sum(sum(gauss));

% return;
end
