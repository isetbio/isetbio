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
   % family, size, dpi
   font = fontCreate('g', 'Georgia', 24, 96); 
   vcNewGraphWin; imagesc(font.bitmap);
%}

if notDefined('letter'), letter = 'g'; end
if notDefined('sz'),     sz = 14; end
if notDefined('family'), family = 'Georgia'; end
if notDefined('dpi'),    dpi = 96; end

font.type = 'font';
font.name = lower(sprintf('%s-%s-%i-%i', letter, family, sz, dpi));
font.character = letter;
font.size = sz;
font.family = family;
font.style = 'NORMAL';
font.dpi = dpi;

% Reads the cached font bit maps.  Microsoft had proprietary technology for
% generating fonts, and we stored a number of example.
font.bitmap = fontBitmapGet(font);
    
end

