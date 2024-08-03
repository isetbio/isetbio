function [varargout] = CalibrateFitLinMod(calOrCalStruct)
% Usage (new style):
% CalibrateFitLinMod(calOrCalStruct)
% where calOrCalStruct is a @CalStruct object.
%
% Usage (old style):
% cal = CalibrateFitLinMod(calOrCalStruct)
% where cal and calOrCalStruct are both old-style cal structures.
% 
%
% Fit the linear model to spectral calibration data.
%
% 3/26/02   dhb  Pulled out of CalibrateMonDrvr.
% 3/27/02   dhb  Add case of nPrimaryBases == 0.
% 2/15/10   dhb  Fix so that first basis vector is good approximation to max
%               input primary spectrum.
%           dhb  Normalize basis vectors so that their max power matches that 
%               of first component.
% 4/30/10   dhb  Execute yoked fit if yokedGamma flag is set.
% 5/25/10   dhb, ar Change yoked field names to match
% 5/26/10   dhb, ar Still fussing with names.
% 5/28/10   dhb, ar Pull out yoked fitting from here -- too confusing.
% 5/27/12   dhb     Handle case where there are more measurements than wavelength samples.
% 5/28/14   npc  Modifications for accessing calibration data using a @CalStruct object.
%                The first input argument can be either a @CalStruct object (new style), or a cal structure (old style).
%                Passing a @CalStruct object is the preferred way because it results in 
%                (a) less overhead (@CalStruct objects are passed by reference, not by value), and
%                (b) better control over how the calibration data are accessed.


% Specify @CalStruct object that will handle all access to the calibration data.
[calStructOBJ, inputArgIsACalStructOBJ] = ObjectToHandleCalOrCalStruct(calOrCalStruct);
if (inputArgIsACalStructOBJ)
    % The input (calOrCalStruct) is a @CalStruct object. Make sure that CalibrateFitLinMod is called with no return variables.
    if (nargout > 0)
        error('There should be NO return parameters when calling CalibrateFitLinMod with a @CalStruct input. For more info: doc CalibrateFitLinMod.');
    end
else
    % The input (calOrCalStruct) is a cal struct. Clear it to avoid  confusion.
    clear 'calOrCalStruct';
end
% From this point onward, all access to the calibration data is accomplised via the calStructOBJ.

% Extract needed data
S               = calStructOBJ.get('S');
nDevices        = calStructOBJ.get('nDevices');
nPrimaryBases   = calStructOBJ.get('nPrimaryBases');
nMeas           = calStructOBJ.get('nMeas');
mon             = calStructOBJ.get('mon');

Pmon = zeros(S(3), nDevices*nPrimaryBases);
mGammaRaw = zeros( nMeas, nDevices* nPrimaryBases);
monSVs = zeros(min([ nMeas S(3)]), nDevices);

for i = 1:nDevices
    tempMon = reshape(mon(:,i), S(3), nMeas);
    monSVs(:,i) = svd(tempMon);

    % Build a linear model
    if (nPrimaryBases ~= 0)
        % Get full linear model
        %
        % If we wanted to add the option of ignoring the lower input level
        % measurements when creating the linear model, we think we'd add a
        % parameter and then modify this to:
        %            tempMon(:,nIgnoreForLinMod:end)
        %            monW = monB\tempMon;
        % Plus a similar change as indicated below when the residual basis
        % functions are computed.  But using nPrimaryBases == 0 is
        % basically as good or better than doing this, so we're leaving it
        % alone for now.
        [monB,monW] = FindLinMod(tempMon, nPrimaryBases);

        % Express max measurement within the full linear model.
        % This is the first basis function.
        tempB = monB*monW(:, nMeas);
        maxPow = max(abs(tempB));

        % Get residuals with respect to first component
        residMon = tempMon-tempB*(tempB\tempMon);

        % If linear model dimension is greater than 1,
        % fit linear model of dimension-1 to the residuals.
        % Take this as the higher order terms of the linear model.
        %
        % Also normalize each basis vector to max power of first
        % component, and make sure convention is that this max
        % is positive.
        if (nPrimaryBases > 1)
            % Change to FindLinMod(residMon(:,nIgnoreForLinMod:end, nPrimaryBases-1)
            residB = FindLinMod(residMon, nPrimaryBases-1);
            for j = 1: nPrimaryBases-1
                residB(:,j) = maxPow*residB(:,j)/max(abs(residB(:,j)));
                [nil,index] = max(abs(residB(:,j)));
                if (residB(index,j) < 0)
                    residB(:,j) = -residB(:,j);
                end
            end
            monB = [tempB residB];
        else
            monB = tempB;
        end

        % Zero means build one dimensional linear model just taking max measurement
        % as the spectrum.
    else
        calStructOBJ.set('nPrimaryBases', 1);
        monB = tempMon(:, nMeas);
    end

    % Find weights with respect to adjusted linear model and
    % store
    monW = FindModelWeights(tempMon,monB);
    for j = 1:nPrimaryBases
        mGammaRaw(:,i+(j-1) * nDevices) = (monW(j,:))';
        Pmon(:,i+(j-1) * nDevices) = monB(:,j);
    end
end

% Update calibration structure.               
calStructOBJ.set('P_device', Pmon);              
calStructOBJ.set('T_device', WlsToT(S));          
calStructOBJ.set('rawGammaTable', mGammaRaw);     
calStructOBJ.set('monSVs', monSVs);               

if (~inputArgIsACalStructOBJ)
    % Old-style functionality. Return modified cal.
    varargout{1} = calStructOBJ.cal;
    % calStructOBJ is not needed anymore. So clear it from the memory.
    clear 'calStructOBJ'
end
