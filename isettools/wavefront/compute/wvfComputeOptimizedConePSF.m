function wvfParams = wvfComputeOptimizedConePSF(wvfParams)
% wvfParams = wvfComputeOptimizedConePSF(wvfParams)
%
% Optimize the PSF seen by the cones, given the cone sensitivities, a
% weighting spectral power distribution, and a criterion.  Optimization is
% performed on the defocus parameter
%
% 8/26/11  dhb  Wrote it.
% 8/29/11  dhb  Don't need to center or circularly average here.
%          dhb  Print warning if optimal value is at search bound.
% 9/7/11   dhb  Rename.  Use wvfParams for i/o.

options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','active-set');
%options = optimset(options,'TypicalX',0.1,'DiffMinChange',1e-3);

% If the parallel toolbox is present, check and use if available
if exist('matlabpool','builtin')
    if (exist('IsCluster','file') && IsCluster && matlabpool('size') > 1)
        options = optimset(options,'UseParallel','always');
    end
end

% Initial defocus and bounds (diopters)
diopterBound = 4;
x0 = 0;
vlb = -diopterBound;
vub = -vlb;

% Optimize focus
x = fmincon(@InlineMinFunction,x0,[],[],[],[],vlb,vub,[],options);

% Set up return values
defocusDiopters = x;
if (abs(defocusDiopters) >= diopterBound)
    fprintf('WARNING: defocus found in wvfComputeOptimizedConePSF is at search limit of %0.1f diopters\n',diopterBound)
end
[f,tmpWvfParams] = InlineMinFunction(defocusDiopters);
wvfParams = tmpWvfParams;

    function [f,tmpWvfParams] = InlineMinFunction(x)
        tmpWvfParams = wvfParams;
        tmpWvfParams.defocusDiopters = x;
        tmpWvfParams = wvfComputeConePSF(tmpWvfParams);
        nCones = size(tmpWvfParams.T_cones,1);
        f = 0;
        for j = 1:nCones
            %temppsf = psfCircularlyAverage(psfCenter(conepsf(:,:,j)));
            critRadius(j) = psfFindCriterionRadius(tmpWvfParams.conepsf(:,:,j),tmpWvfParams.criterionFraction);
            f = f + tmpWvfParams.coneWeights(j)*critRadius(j);
        end
        
    end
end
