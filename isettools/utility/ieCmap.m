function cm = ieCmap(cName, num, gam)
% Prepare simple color maps
%
% Syntax:
%   cm = ieCmap(cName, [num], [gam])
%
% Description:
%    Create a simple color map using the color map type, the number of
%    colors, and the luminance.
%
% Inputs:
%    cName - (Optional) Color map name. Default 'rg'. Options are:
%               {'rg', 'redgreen'}
%               {'by', 'blueyellow'}
%               {'bw', 'blackwhite', 'luminance'}
%    num   - (Optional) Number of elements in the color map. Default 256.
%    gam   - (Optional) Set gamma value (default = 1)
%
% Outputs:
%    cm    - The color map to return
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/10       Copyright ImagEval Consultants, LLC, 2010
%    11/30/17  jnm  Formatting
%    01/18/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    rg  = ieCmap('rg', 256,1);
    imagesc(rand(128,128)); colormap(rg); colorbar;
    rg  = ieCmap('rg', 256, .5); colormap(rg); colorbar;
%}
%{
	by  = ieCmap('by', 256);
    plot(by)
    imagesc(rand(128,128)); colormap(by); colorbar;
%}
%{
	lum = ieCmap('bw', 256, 0.3);
    plot(lum)
%}

if notDefined('cName'), cName = 'rg'; end
if notDefined('num'), num = 256; end
if notDefined('gam'), gam = 1; end

cName = ieParamFormat(cName);

switch cName
    case {'redgreen', 'rg'}
        a = linspace(0, 1, num);
        cm = [a(:), flipud(a(:)), 0.5 * ones(size(a(:)))];
        
    case {'blueyellow', 'by'}
        a = linspace(0, 1, num);
        cm = [a(:), a(:), flipud(a(:))];
        
    case {'luminance', 'blackwhite', 'bw'}
        cm = gray(num);
        
    otherwise
        error('Unknown color map name %s\n', cName);
end

% All maps can have the gamma applied
cm = cm .^ gam;

end
