function tickLocs = ieChooseTickMarks(val, nTicks)
% Choose sensible values for tick marks (or at least try)
%
% Syntax:
%   tickLocs = ieChooseTickMarks(val, [nTicks])
%
% Description:
%    Choose tick mark locations on an axis with the values in val.
%
%    There are examples in the code. To access the examples, type 'edit
%    ieChooseTickMarks.m' into the Command Window.
%
%    The code below contains examples of function usage. To access, type
%    'edit ieChooseTickMarks.m' into the Command Window.
%
% Inputs:
%    val      - Integer Array. Axis values.
%    nTicks   - (Optional) Integer. Number of ticks. Default = 10.
%
% Outputs:
%    tickLocs - The locations of the tick marks.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * This routine fails when nTicks > (mx - mn)
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    02/28/18  jnm  Formatting


% Examples:
%{
    tickLocs = ieChooseTickMarks(0:500, 5)
%}

if notDefined('nTicks'), nTicks = 10; end

% Choose a reasonable number of positive tick marks
% [Note: JNM - Contrary to the statement above, these values are not
% guaranteed to be positive.]
mx = round(max(val));
mn = round(min(val));
tickSpacing = round((mx - mn) / nTicks);

E = floor(log10(mx - mn));
S = 10 ^ (E - 1);
tickSpacing = S * round(tickSpacing / S);

if tickSpacing == 0
    error('Tick spacing cannot be zero.')
end

if mn < 0 && mx > 0
    tickLocs = 0:tickSpacing:mx;
    tickLocs = [fliplr(-1 * (tickSpacing:tickSpacing:abs(mn))), tickLocs];
else
    tickLocs = mn:tickSpacing:mx;
end

end