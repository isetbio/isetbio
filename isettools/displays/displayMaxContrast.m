function alpha = displayMaxContrast(signalDir, backDir)
% Find scalar that produces max contrast for signal in given background
%
% Syntax:
%   alpha = displayMaxContrast(signalDir, backDir)
%
% Description:
%    Find the scalar the produces the maximum contrast for a signal
%    given a background.
%
%    With the scalar alpha applied to the signal, the total signal
%    reaches one or the other display boundary in the equation.
%    This routine is useful, for example, in finding the maximum
%    contrast cone-isolating stimulus that can be displayed on a
%    particular monitor.
%
%        0 <= alpha * signalDir + backDir
%             alpha * signalDir + backDir <= 1
%
% Inputs:
%    signalDir - The signal in question.
%    backDir   - Provided background
%
% Outputs:
%    alpha     - Numeric. The desired scalar.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    05/16/18  jnm  Formatting

for ii = 1:3
    if signalDir(ii) > 0
        mx(ii) = (1 - backDir(ii)) / signalDir(ii);
    else
        mx(ii) = abs(-backDir(ii) / signalDir(ii));
    end
end
alpha = min(mx);

return
