function font = fontCreate(letter, family, sz, dpi)
% Create a font structure
%
% Syntax:
%   font = fontCreate([letter], [family], [sz], [dpi])
%
% Description:
%    Create a font structure. There are a number of parameters set as
%    defaults, in addition to the options you may specify when creating the
%    font.
%
%    For example, the bitmap of the font is set so the letter is black on a
%    white background. Using the setFont function you may use 1 - bitmap to
%    flip it to white on a black background.
%
%    There are examples contained in the code below. To access the
%    examples, type 'edit fontCreate.m' into the Command Window.
%
% Inputs:
%    letter - (Optional) Char. The letter. Default 'g'.
%    sz     - (Optional) Integer. Font size. Default 14.
%    family - (Optional) String. Font family. Default 'Georgia'
%    dpi    - (Optional) Integer. Dots per inch. Default 96.
%
% Outputs:
%    font   - The created font structure.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * TODO: create way to read cached fonts (See note below).

% History:
%    xx/xx/14  BW   Vistasoft group, 2014
%    02/27/18  jnm  Formatting

% Examples:
%{
   font = fontCreate;
   font = fontCreate('A', 'Georgia', 24, 96); 
   vcNewGraphWin;
    imagesc(font.bitmap);

   font = fontCreate('l', 'georgia', 14, 72);
%}

if notDefined('letter'), letter = 'g'; end
if notDefined('sz'), sz = 14; end
if notDefined('family'), family = 'Georgia'; end
if notDefined('dpi'), dpi = 96; end

font.type = 'font';
font.name = lower(sprintf('%s-%s-%i-%i', letter, family, sz, dpi));
font.character = letter;
font.size = sz;
font.family = family;
font.style = 'NORMAL';
font.dpi = dpi;

% [Note: XXX - Need to make a way to read the cached fonts rather than the
% method here.]
% The dpi will be added and we will read the overscaled (x) fonts. Then we
% will put them in the display structure, and maybe filter if we decide to.
font.bitmap = fontBitmapGet(font);
    
end

%%
function bitmap = fontBitmapGet(font)
% Read a bitmap from a stored file (data/fonts) or create one on the screen
%
% Syntax:
%   bitmap = fontBitmapGet(font)
%
% Description:
%    This function helps get the font bitmap from the system. The basic
%    idea of this function is to draw a letter on a canvas in the
%    background and read back the frame
%
% Inputs:
%   font    - font data structure (see fontCreate)
%
% Outputs:
%    bitMap - bitmap image
%
% Optional key/value pairs:
%    None.
%

% History:
%    05/xx/14  HJ   May, 2014
%    02/27/18  jnm  Formatting

% Examples:
%{
   font = fontCreate;
    bitmap = fontBitmapGet(font);
%}

%% Init
if notDefined('font'), error('font required'); end

try
    % See if we have a font stored for this variable
    name = fontGet(font, 'name');
    fName = sprintf('%s.mat', name);
    fName = fullfile(isetbioDataPath, 'fonts', fName);
    load(fName);
    b = bmSrc.dataIndex;     
    b = 1 - b;  %Make black on white
    padsize = 3 * ceil(size(b, 2) / 3) - size(b, 2);
    b = padarray(b, [0 padsize], 1, 'pre');

    % Reformat the bit maps, which are not always a multiple of 3 wide
    bitmap = ones(size(b, 1), ceil(size(b, 2) / 3), 3);
    for ii = 1:3, bitmap(:, :, ii) = b(:, ii:3:end); end

    return;
catch
    
    % No font file found. Create an example on the screen
    character = fontGet(font, 'character');
    fontSz = fontGet(font, 'size');
    family = fontGet(font, 'family');
    
    %% Set up canvas
    hFig = figure('Position', [50 50 150+fontSz 150+fontSz], ...
        'Units', 'pixels', ...
        'Color', [1 1 1], ...
        'Visible', 'off');
    axes('Position', [0 0 1 1], 'Units', 'Normalized');
    axis off;

    %% Draw character and capture the frame
    % Place each character in the middle of the figure
    texthandle = text(0.5, 1, character, ...
        'Units', 'Normalized', ...
        'FontName', family, ...
        'FontUnits', 'pixels', ...
        'FontSize', fontSz, ...
        'HorizontalAlignment', 'Center', ...
        'VerticalAlignment', 'Top', ...
        'Interpreter', 'None', ...
        'Color', [0 0 0]);
    drawnow;

    % Take a snapshot
    try
        bitMap = hardcopy(hFig, '-dzbuffer', '-r0');
    catch
        bitMap = getframe(hFig);
        bitMap = bitMap.cdata;
    end

    delete(texthandle);

    % Crop height as appropriate
    bwBitMap = min(bitMap, [], 3);
    bitMap = bitMap(find(min(bwBitMap, [], 2) < 255, 1, 'first') : ...
        find(min(bwBitMap, [], 2) < 255, 1, 'last'), :, :);

    % Crop width to remove all white space
    bitMap = bitMap(:, find(min(bwBitMap, [], 1) < 255, 1, 'first') : ...
        find(min(bwBitMap, [], 1) < 255, 1, 'last'), :);

    % Invert and store in binary format
    bitmap = zeros(size(bitMap));
    bitmap(bitMap > 127) = 1;

    close(hFig);

end

end
