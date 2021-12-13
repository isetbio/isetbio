function XYZ = LjgToXYZ(Ljg)
% XYZ = LjgToXYZ(Ljg)
%
% Convert OSA Ljg to XYZ (10 degree).  Works by using numerical
% search to invert XYZToLjg.  See XYZToLjg for details on
% formulae used.  See also OSAUCSTest.
%
% This can return imaginary values if you pass XYZ values
% that are outside reasonable physical gamut limits.
%
% 3/27/01  dhb      Wrote it.
% 3/4/05   dhb      Handle new version of optimization toolbox, too.
% 9/23/12  dhb, ms  Update options for current Matlab versions.

% On Octave we need the 'optim' package v1.6.0 or later for fmincon() support,
% and Octave 6 or later for support for handles to nested functions like @LjgToXYZFun below:
if IsOctave
    v = version;
    if str2num(v(1)) < 6
        error('For use with Octave, you need at least Octave version 6.');
    end

    try
        % Try loading the optim package with the optimization functions:
        pkg load optim;
    catch
        error('For use with Octave, you must install the ''optim'' package from Octave Forge. See ''help pkg''.');
    end

    % Got optim package loaded. Does it support fmincon()?
    if ~exist('fmincon')
        error('For use with Octave, you need at least version 1.6.0 of the ''optim'' package from Octave Forge.');
    end
end



% Check for needed optimization toolbox, and version.
if (exist('fmincon') == 2)
    % Search options
    options = optimset('fmincon');
    options = optimset(options, 'Display', 'off', 'Algorithm', 'active-set');

    % Search bounds -- XYZ must be positive.
    vlb = zeros(3,1);
    vub = [200 100 200]';

    % Do the search for each passed value.
    XYZ = zeros(size(Ljg));
    n = size(XYZ,2);
    for i = 1:n
        x0 = xyYToXYZ([.28 .28 30]');
        LjgCurrent = Ljg(:,i);
        x = fmincon(@LjgToXYZFun,x0,[],[],[],[],vlb,vub,[],options);
        XYZ(:,i) = x;
    end
elseif (exist('constr') == 2)
    % Search options
    options = foptions;
    options(1) = 0;

    % Search bounds -- XYZ must be positive.
    vlb = zeros(3,1);
    vub = [200 100 200]';

    % Do the search for each passed value.
    XYZ = zeros(size(Ljg));
    n = size(XYZ,2);
    for i = 1:n
        x0 = xyYToXYZ([.28 .28 30]');
        LjgCurrent = Ljg(:,i);
        x = constr(@LjgToXYZFun,x0,options,vlb,vub);
        XYZ(:,i) = x;
    end
else
    error('LjgToXYZ requires the optional Matlab Optimization Toolbox from Mathworks');
end

        function [f,g] = LjgToXYZFun(XYZ)
            % [f,g] = LjgToXYZFun(XYZ)
            %
            % Nested optimization function for LjgToXYZ numerical
            % search. Not vectorized. Uses LjgCurrent from parent function.
            %
            % 3/27/01  dhb  Wrote it.
            % 12/4/20  mk   Turned it into nested function.

            % Convert back.
            Ljg1 = XYZToLjg(XYZ);

            % Handle case where XYZ is so weird that
            % an imaginary value is returned.
            if (any(~isreal(Ljg1)))
                sdiff = (abs(Ljg1)-LjgCurrent).^2;
                f = sum(sdiff);
                g = 10;
            else
                sdiff = (Ljg1-LjgCurrent).^2;
                f = sum(sdiff);
                g = -XYZ;
            end
        end
end
