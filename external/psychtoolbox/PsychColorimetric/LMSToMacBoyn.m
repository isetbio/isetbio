function [ls, factorsLMS] = LMSToMacBoyn(LMS,T_cones,T_lum,lumReturnFlag)
% [ls, factorsLMS] = LMSToMacBoyn(LMS,[T_cones,T_lum],[lumReturnFlag)
%
% Compute MacLeod-Boynton chromaticity from cone exciation coordinates.
% This is L/Lum and S/Lum, with appropriate normalization as described
% below. Recommended usage yields MB chromacities according to CIE 170-2:2015.
%
% ** Recommended Usage: Compute LMS with respect to some specified T_cones,
% and pass both T_cones as well as the cooresponding T_lum (the photopic
% luminosity function.) T_lum should be a linear combination of the L and M
% cone fundamentals.
% 
% This routine will then scale passed L and M values so that they sum to
% the best linear approximation of luminance and then normalize the L
% excitation by luminance (as computed by the linear combination) to obtain
% the l chromaticity. It will normalize the S excitaiton by luminance to
% obtain s chromaticity, with an overall scaling so that the maximum value
% of this chromaticity is 1 taken over the visible spectrum.
%
% Note that the s cone scaling can vary a bit depending on the wavelength
% sampling of the passed T_cones and T_lum, since the max is taken over
% these. If you use the T_cones_ss2/T_cones_ss10 and T_CIE_Y2/T_CIE_Y10
% files provided in PTB, the default sampling is at 1 nm and this is fine.
% If you use subsampled wavelength spacing, the computation of the s cone
% scaling will begin to deviate from the standard.  But so will your
% computation of LMS values, so this isn't really an issue specific to this
% routine.
%
% When you use the CIE cone fundamentals and corresponding luminance
% functions, this procedure yields the MacLeod-Boynton chromaticity
% diagrams as specified in CIE 170-2:2015.
%
% It would have been smarter to write this routine to return lsY rather
% than just ls, long ago.  Getting Y back makes it easier to convert back
% to LMS, and also be consistent with the way other conversions to
% chromaticity are set up.  Changing this now would probably break some
% existing code, so instead this routine now takes an optional fourth
% argument that returns the luminance used in the denominator.  If you ask
% for this, you need to provide T_LMS and T_Y as well, which is recommended
% for clarity in any case.  Set lumReturnFlag = 1 on call to get back three
% vectors with Y (the normalizing denominator) as the third coordinate.  In
% general this will not be the same as L+M in the passed cone system, but
% rather a weighted combination of L and M, with the weights computed from
% the passed T_LMS and T_Y.
%
% The factorsLMS return value is a column vector with three entries that
% specifies how this routine decided to scale the passed LMS cone
% fundamentals to best approximate luminance and get the s axis scaling.
% Note that if you pass LMS as the empty matrix, you can obtain this vector without doing much else.
%
% ** Legacy Usage: Just pass LMS values. In this case, we assume that the
% passed LMS values were computed with respect to the Smith-Pokorny
% fundamentals normalized to a peak of 1 and Judd-Vos luminance (more or
% less).  That is, this usage assumes LMS was computed using the
% fundamentals stored in PTB's T_cones_sp. This is old usage and preserved
% for backwards compatibility, but the three argument usage as described
% above is preferred for clarity. Moreover, in this case, the s
% chromaticity is not further normalized.  This leads to S chromaticities
% considerably larger than those obtained with the new usage.
%
% ** A Backwards Incompatibility. The scaling for s chromaticity to match
% CIE 170-2:2015 was introduced in Janurary 2019 and is not backwards
% compatible with previous behavior when T_cones and T_lum are passed.
% Preserving such compatibility did not seem important enough relative to
% the gains of having this work as now specified in the CIE standard.
%
% 10/30/97  dhb  Wrote it.
% 7/9/02    dhb  T_cones_sp -> T_cones on line 20.  Thanks to Eiji Kimura.
% 1/23/19   dhb  Scale s chromaticity value to be consistent with CIE
%                170-2:2015, when T_cones and T_lum are passed.  This is
%                not backwards compatible with previous scaling, but it
%                seems good to match the standard. Thanks to Danny Garside
%                for pointing out the scaling specified in the 2015
%                standard.
% 06/16/25 dhb   Provide lurReturnFlag and corresponding output option, but
%                keep code backward compatible when that is not passed.
%          dhb   Return the factorsLM column vector computed by this
%                routine, so that it is easily available to the caller if 
%                wanted.

% Examples:
%{
    % Recreate the spectrum locus and equal energy white shown in Figure 8.2
    % of CIE 170-2:2015. Also performs a regression check.
    clear; close all;
    load T_cones_ss2
    load T_CIE_Y2
    lsSpectrumLocus = LMSToMacBoyn(T_cones_ss2,T_cones_ss2,T_CIE_Y2);
 
    % Compute the sum of the ls values in the spectrum locus, and compare
    % to the value that this example computed in February 2019, entered
    % here to four places as 412.2608.  This comparison provides a
    % check that this routine still works the way it did when we put in the
    % check.
    check = round(sum(lsSpectrumLocus(:)),4);
    if (abs(check-412.2608) > 1e-4)
        error('No longer get same check value as we used to');
    end

    % Compute ls for equal energy white
    LMSEEWhite = sum(T_cones_ss2,2);
    lsEEWhite = LMSToMacBoyn(LMSEEWhite,T_cones_ss2,T_CIE_Y2); 

    % Plot
    figure; hold on;
    plot(lsSpectrumLocus(1,:)',lsSpectrumLocus(2,:)','r','LineWidth',3);
    plot(lsEEWhite(1),lsEEWhite(2),'bs','MarkerFaceColor','b','MarkerSize',12);
    xlim([0.4 1]); ylim([0,1]);
    xlabel('l chromaticity'); ylabel('s chromaticity');
    title('CIE 170-2:2015 Figure 8.2');
%}
%{
    % Return the factorsLMS vector for the caller.
    % This does CIE-land version.
    clear; close all;
    load T_cones_ss2
    load T_CIE_Y2
    [~,factorsLMS] = LMSToMacBoyn([],T_cones_ss2,T_CIE_Y2);
    factorsLMS

    % Also ought to work with MacBoynToLMS but you don't need
    % the flag since this is a new routine and only works one way.
    [~,factorsLMS] = MacBoynToLMS([],T_cones_ss2,T_CIE_Y2);
    factorsLMS

    % This version does the Smith-Pokorny to Judd-Vos luminance scaling
    % explicitly. Gets a different answer in the fourth place than the
    % default version, which is probably due to choice of wavelength
    % sampling or other small numerical difference.
    clear; close all;
    load T_cones_sp
    load T_xyzJuddVos
    [~,factorsLMS] = LMSToMacBoyn([],T_cones_sp,T_xyzJuddVos(2,:));
    factorsLMS

    % This version (no longer recommended) returns the default scaling of 
    % Smith-Pokorny fundamentals to get Judd-Vos luminance
    [~,factorsLMS] = LMSToMacBoyn([]);
    factorsLMS
%}
%{
    % Recreate the diagram with respect to the Smith-Pokorny fundamentals.
    clear; close all;
    load T_cones_sp
    load T_xyzJuddVos
    T_xyzJuddVos = SplineCmf(S_xyzJuddVos,T_xyzJuddVos,S_cones_sp);
    lsSpectrumLocus = LMSToMacBoyn(T_cones_sp,T_cones_sp,T_xyzJuddVos(2,:));

    % Compute ls for equal energy white
    LMSEEWhite = sum(T_cones_sp,2);
    lsEEWhite = LMSToMacBoyn(LMSEEWhite,T_cones_sp,T_xyzJuddVos(2,:)); 

    % Plot
    figure; hold on;
    plot(lsSpectrumLocus(1,:)',lsSpectrumLocus(2,:)','r','LineWidth',3);
    plot(lsEEWhite(1),lsEEWhite(2),'bs','MarkerFaceColor','b','MarkerSize',12);
    xlim([0.4 1]); ylim([0,1]);
    xlabel('l chromaticity'); ylabel('s chromaticity');
    title('MacBoyn wrt Smith-Pokorny/JuddVos');
%}
%{
    % Demonstrate invariance of ls after scaling of cones and luminance, as
    % long as LMS valued are computed with respect to passed cones and
    % luminance.
    clear; close all;
    load T_cones_ss2
    load T_CIE_Y2

    % Spectrum locus ls chromaticity with no scaling.
    lsSpectrumLocus = LMSToMacBoyn(T_cones_ss2,T_cones_ss2,T_CIE_Y2); 

    % Do some scaling and recompute spectrum locus.
    % The choices of 1.8, 3, 0.05 as scaling are
    % just three numbers I made up. You can use any three numbers and it
    % should still work, modulo any numerical issues with very big or
    % very small numbers.
    T_CIE_Y2_scaled = 1.8*T_CIE_Y2;
    T_cones_ss2_scaled = T_cones_ss2;
    T_cones_ss2_scaled(1,:) = 3*T_cones_ss2(1,:);
    T_cones_ss2_scaled(3,:) = 0.05*T_cones_ss2(3,:);
    lsSpectrumLocusScaled = LMSToMacBoyn(T_cones_ss2_scaled,T_cones_ss2_scaled,T_CIE_Y2_scaled); 

    % Make sure the difference is very small numerically.
    diff = max(abs(lsSpectrumLocus(:)-lsSpectrumLocusScaled(:)));
    fprintf('Difference in spectrum locus chromaticity after scaling is %0.5f (should be small)\n',diff);

    % Plot the locus computed both ways to compare visually.
    figure; hold on;
    plot(lsSpectrumLocus(1,:)',lsSpectrumLocus(2,:)','r','LineWidth',3);
    plot(lsSpectrumLocusScaled(1,:)',lsSpectrumLocusScaled(2,:)','g','LineWidth',2);
    xlim([0.4 1]); ylim([0,1]);
    xlabel('l chromaticity'); ylabel('s chromaticity');
%}
%{
    % Recreate the spectrum locus and equal energy white shown in Figure 8.2
    % of CIE 170-2:2015. Also performs a regression check.
    % 
    % This version also returns luminance, which makes it easy to invert
    clear; close all;
    load T_cones_ss2
    load T_CIE_Y2
    lsYSpectrumLocus = LMSToMacBoyn(T_cones_ss2,T_cones_ss2,T_CIE_Y2,1);
 
    % Compute the sum of the ls values in the spectrum locus, and compare
    % to the value that this example computed in February 2019, entered
    % here to four places as 412.2608.  This comparison provides a
    % check that this routine still works the way it did when we put in the
    % check.
    temp = lsYSpectrumLocus(1:2,:);
    check = round(sum(temp(:)),4);
    if (abs(check-412.2608) > 1e-4)
        error('No longer get same check value as we used to');
    end

    % Compute ls for equal energy white
    LMSEEWhite = sum(T_cones_ss2,2);
    lsYEEWhite = LMSToMacBoyn(LMSEEWhite,T_cones_ss2,T_CIE_Y2,1); 

    % Let's go back to LMS
    LMSSpectrumLocusCheck = MacBoynToLMS(lsYSpectrumLocus,T_cones_ss2,T_CIE_Y2);
    LMSEEWhiteCheck = MacBoynToLMS(lsYEEWhite,T_cones_ss2,T_CIE_Y2);

    % Check self inversion
    if (max(abs(LMSSpectrumLocusCheck(:) - T_cones_ss2(:))) > 1e-6)
        error('Self inversion failure on spectrum locus');
    end
    if (max(abs(LMSEEWhiteCheck(:) - LMSEEWhite(:))) > 1e-4)
        error('Self inversion failure on spectrum locus');
    end 
%}

% Scale LMS so that L+M = luminance and S cone value corresponds to a
% fundamental with a max of 1.
validInput = false;
if (nargin == 1)
    factorsLM = [0.6373 0.3924]';
    factorS = 1;
    lumReturnFlag = 0;
    validInput = true;
elseif (nargin >= 3 )
	factorsLM = (T_cones(1:2,:)'\T_lum');
    factorS = 1/max(T_cones(3,:)./(factorsLM(1)*T_cones(1,:) + factorsLM(2)*T_cones(2,:)));
    validInput = true;
end
factorsLMS = [factorsLM ; factorS];

if (nargin == 3)
    lumReturnFlag = 0;
elseif (nargin > 4)
    validInput = false;
end
if (~validInput)
    error('Number of input arguments should be either 1, 3, or 4');
end

% Only do the rest if a non-empty LMS was passed.
% Otherwise return an emtpy matrix as ls
if (~isempty(LMS))
    % Set up return value
    if (lumReturnFlag == 0)
        outputDimension = 2;
    elseif (lumReturnFlag == 1)
        outputDimension = 3;
    else
        error('Bad value of lumReturnFlag passed');
    end

    % Scale LMS so that L+M is our approximation of luminance
    LMS = diag(factorsLMS)*LMS;

    % Compute ls coordinates from LMS
    n = size(LMS,2);
    ls = zeros(outputDimension,n);
    Y = [1 1 0]*LMS;
    ls(1:2,:) = LMS([1 3],:) ./ ([1 1]'*Y);

    % IF we are returning luminance, do so
    if (lumReturnFlag == 1)
        ls(3,:) = Y;
    end
else
    % Handle empty input LMS gracefully
    ls = [];
end

