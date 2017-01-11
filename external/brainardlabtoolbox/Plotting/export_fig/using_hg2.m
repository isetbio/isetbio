%USING_HG2 Determine if the HG2 graphics pipeline is used
%
%   tf = using_hg2(fig)
%
%IN:
%   fig - handle to the figure in question.
%
%OUT:
%   tf - boolean indicating whether the HG2 graphics pipeline is being used
%        (true) or not (false).
%
% 12/19/15  dhb  Modified for 2016b, can't use graphicsversion without warning
function tf = using_hg2(fig)

if (verLessThan('matlab','9.0.0'))
    try
        tf = ~graphicsversion(fig, 'handlegraphics');
    catch
        tf = false;
    end
else
    tf = true;
end
end
