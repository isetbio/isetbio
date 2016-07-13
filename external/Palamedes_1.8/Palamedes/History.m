%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This document will only track changes to the toolbox proper (i.e., those 
%files residing in the 'Palamedes' folder). Changes to files in the 
%PalamedesDemos folder will not be documented here (or elsewhere).
%
%Palamedes: Matlab routines for analyzing psychophysical data.
%
%Nick Prins & Fred Kingdom. palamedes@palamedestoolbox.org
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Version 1.0.0 launch: September 13, 2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Version 1.0.1 release: September 17, 2009
%Modifications:
%
% PAL_SDT_MAFCmatchSample_DiffMod_DPtoPC:
% Modified: Palamedes version 1.0.1 (FK). Changed default value of numReps 
%   from 500000 to 100000
%
% PAL_SDT_MAFCoddity_DPtoPC:
% Modified: Palamedes version 1.0.1 (FK). Changed default value of numReps 
%   from 500000 to 100000
%
% PAL_PFML_BootstrapNonParametricMultiple:
% Modified: Palamedes version 1.0.1 (NP). A warning and suggestion will be 
%   issued when OutOfNum contains ones. 
%
% PAL_PFML_BootstrapNonParametric:
% Modified: Palamedes version 1.0.1 (NP). A suggestion to consider a
%   parametric bootstrap is added to the warning issued when OutOfNum 
%   contains ones.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Version 1.0.2 release: October 1, 2009
%Modifications:
%
% PAL_PFML_GoodnessOfFitMultiple:
% Modified: Palamedes version 1.0.2 (NP). No longer produces 'DivideByZero'
%   warning when B is set to 0 (to avoid the running of simulations).
%
% PAL_PFML_GoodnessOfFit:
% Modified: Palamedes version 1.0.2 (NP). No longer produces 'DivideByZero'
%   warning when B is set to 0 (to avoid the running of simulations).
%
% PAL_PFLR_ModelComparison:
% Modified: Palamedes version 1.0.2 (NP). No longer produces 'DivideByZero'
%   warning when B is set to 0 (to avoid the running of simulations).
% Modified: Palamedes version 1.0.2 (NP). Fixed error in comments section
%   regarding the names of PF routines.
%
% PAL_PFML_Fit:
% Modified: Palamedes version 1.0.2 (NP). Fixed error in comments section
%   regarding the names of PF routines.
%
% PAL_PFML_FitMultiple:
% Modified: Palamedes version 1.0.2 (NP). Fixed error in comments section
%   regarding the names of PF routines.
%
% PAL_PFML_BootstrapNonParametric:
% Modified: Palamedes version 1.0.2 (NP). Fixed error in comments section
%   regarding the names of PF routines.
%
% PAL_PFML_BootstrapParametric:
% Modified: Palamedes version 1.0.2 (NP). Fixed error in comments section
%   regarding the names of PF routines.
%
% PAL_PFML_BootstrapNonParametricMultiple:
% Modified: Palamedes version 1.0.2 (NP). Fixed error in comments section
%   regarding the names of PF routines.
%
% PAL_PFML_BootstrapParametricMultiple:
% Modified: Palamedes version 1.0.2 (NP). Fixed error in comments section
%   regarding the names of PF routines.
%
% PAL_PFML_GoodnessOfFit:
% Modified: Palamedes version 1.0.2 (NP). Fixed error in comments section
%   regarding the names of PF routines.
%
% PAL_PFML_GoodnessOfFitMultiple:
% Modified: Palamedes version 1.0.2 (NP). Fixed error in comments section
%   regarding the names of PF routines.
%
% PAL_Logistic:
% Modified: Palamedes version 1.0.2 (NP). Added some help comments.
%
% PAL_Gumbel:
% Modified: Palamedes version 1.0.2 (NP). Added some help comments.
%
% PAL_HyperbolicSecant:
% Modified: Palamedes version 1.0.2 (NP). Added some help comments.
%
% PAL_CumulativeNormal:
% Modified: Palamedes version 1.0.2 (NP). Added some help comments.
%
% PAL_Weibull:
% Modified: Palamedes version 1.0.2 (NP). Added some help comments.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Version 1.1.0 release: November 23, 2009
%
%Summary of major change: Version 1.1.0 adds the option to custom-define 
%   constraints on the parameters of PFs between several data sets while 
%   simultaneously fitting PFs to multiple data sets. This option is also 
%   added to model comparison, multi-condition bootstrap, and 
%   multi-condition goodness-of-fit routines (specifically: 
%   PAL_PFML_FitMultiple, PAL_PFML_BootstrapNonParametricMultiple,
%   PAL_PFML_BootstrapParametricMultiple, PAL_PFML_GoodnessOfFitMultiple,
%   PAL_PFLR_ModelComparison). All previous functionality of these 
%   routines is retained. Modified and added functions are listed below.
%
%Version 1.1.0 incorporates some other changes also, but these are minor. 
%   All changes are listed below.
%
%Modifications:
%
% PAL_PFML_rangeTries:
% Modified: Palamedes version 1.1.0 (NP): Modified to assign zeros to 
%   entries in multiplier array which correspond to custom-parametrized 
%   parameters.
%
% PAL_PFML_FitMultiple:
% Modified: Palamedes version 1.1.0 (NP). Modified to allow custom-defined
%   reparametrizations of parameters. Also returns the number of free
%   parameters.
%
% PAL_PFML_TtoP:
% Modified: Palamedes version 1.1.0 (NP): Modified to accept custom-defined
%   reparametrizations also.
%
% PAL_PFML_PtoT:
% Modified: Palamedes version 1.1.0 (NP): Modified to accept custom-defined
%   reparametrizations also.
%
% PAL_Entropy:
% Modified: Palamedes version 1.1.0 (NP): upon completion returns all 
%   warning states to prior settings.
%
% PAL_MLDS_Bootstrap:
% Modified: Palamedes version 1.1.0 (NP): upon completion returns all 
%   warning states to prior settings.
%
% PAL_PFLR_ModelComparison:
% Modified: Palamedes version 1.1.0 (NP). Modified to allow custom-defined
%   reparametrization of parameters.
% Modified: Palamedes version 1.1.0 (NP): upon completion returns all 
%   warning states to prior settings.
%
% PAL_PFLR_TLR:
% Modified: Palamedes version 1.1.0 (NP). Modified to accept
%   custom-reparametrization of parameters.
%
% PAL_PFML_BootstrapNonParametric:
% Modified: Palamedes version 1.1.0 (NP): upon completion returns all 
%   warning states to prior settings.
%
% PAL_PFML_BootstrapNonParametricMultiple:
% Modified: Palamedes version 1.1.0 (NP). Modified to allow custom-defined
%   reparametrization of parameters.
% Modified: Palamedes version 1.1.0 (NP): upon completion returns all 
%   warning states to prior settings.
%
% PAL_PFML_BootstrapParametric:
% Modified: Palamedes version 1.1.0 (NP): upon completion returns all 
%   warning states to prior settings.
%
% PAL_PFML_BootstrapParametricMultiple:
% Modified: Palamedes version 1.1.0 (NP). Modified to allow custom-defined
%   reparametrization of parameters.
% Modified: Palamedes version 1.1.0 (NP): upon completion returns all 
%   warning states to prior settings.
%
% PAL_PFML_GoodnessOfFit:
% Modified: Palamedes version 1.1.0 (NP): upon completion returns all 
%   warning states to prior settings.
%
% PAL_PFML_GoodnessOfFitMultiple:
% Modified: Palamedes version 1.1.0 (NP): upon completion returns all 
%   warning states to prior settings.
% Modified: Palamedes version 1.1.0 (NP). Modified to accept custom-defined
%   reparametrizations of parameters.
%
% PAL_PFML_LLNonParametric:
% Modified: Palamedes version 1.1.0 (NP). Returns the number of free
%   parameters.
%
%Added Routines:
%
% PAL_whatIs:
% Introduced: Palamedes version 1.1.0 (NP): Determines variable type. 
%   Internal function.
%
% PAL_PFML_IndependentFit:
% Introduced: Palamedes version 1.1.0 (NP): Determines whether PFs to 
%   multiple conditions can be fit individually or whether 
%   interdependencies exist. Internal function.
%
% PAL_PFML_LLsaturated:
% Introduced: Palamedes version 1.1.0 (NP): Returns Log Likelhood and 
%   number of parameters in saturated model.
%
% PAL_PFML_setupParametrizationStruct:
% Introduced: Palamedes version 1.1.0 (NP): Creates a parameter 
%   reparametrization structure for (optional) use in functions which allow
%   specification of a model regarding the parameters of PFs across several
%   datasets.
%
% PAL_PFML_CustomDefine:
% Introduced: Palamedes version 1.1.0 (NP): This file only contains 
%   instructions on the use of custom reparametrization.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Version 1.1.1 Release October 31, 2010
%
% 
%
% Changed comments in PAL_AMPM_Demo to make explicit that values for Slope
%   estimate are in log (base 10) units.
%
% PAL_AMPM_setupPM:
% Modified: Palamedes version 1.1.1 (NP): Included previously omitted 
%   lines:  
%   PM.numTrials = 50;
%   PM.response = [];
%
%   Made explicit in the help section that log values used are log base 10.
%
% PAL_Weibull:
% Modified: Palamedes version 1.1.1 (NP). Allowed gamma and lambda to be
%   multidimensional arrays.
%
% PAL_HyperbolicSecant:
% Modified: Palamedes version 1.1.1 (NP). Allowed gamma and lambda to be
%   multidimensional arrays.
%
% PAL_CumulativeNormal:
% Modified: Palamedes version 1.1.1 (NP). Allowed gamma and lambda to be
%   multidimensional arrays.
%
% PAL_Logistic:
% Modified: Palamedes version 1.1.1 (NP). Allowed gamma and lambda to be
%   multidimensional arrays.
%
% Added routines:
%
% PAL_findMax:
% Introduced: Palamedes version 1.1.1 (NP): find value and position of 
%   maximum in 2 or 3D array
%
% PAL_PFML_paramsTry:
% Introduced: Palamedes version 1.1.1 (NP): Generate jitter on values of 
%   guesses to be supplied to PAL_PFML_Fit or PAL_PFML_FitMultiple as 
%   initial values in search.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Version 1.2.0 Release March 23, 2011
%
% 
%
% Changed PAL_PFML_Demo to include new features detailed below
% Introduced PAL_PFML_BruteForceInitials_Demo to demonstrate new features
% detailed below.
%
% PAL_AMPM_CreateLUT:
% Modified: Palamedes version 1.2.0 (NP). Corrected error in function name.
%
% PAL_AMPM_setupPM:
% Modified: Palamedes version 1.2.0 (NP): Fixed some stylistic nasties
%   (failing to pre-allocate arrays, etc.)
%
% PAL_AMRF_setupRF:
% Modified: Palamedes version 1.2.0 (NP): Fixed some stylistic nasties
%   (failing to pre-allocate arrays, etc.)
%
% PAL_AMUD_setupUD:
% Modified: Palamedes version 1.2.0 (NP) Added Garia-Perez reference in
%   comments.
%
% PAL_CumulativeNormal:
% Modified: Palamedes version 1.2.0 (NP). Added inverse PF and derivative of
%   PF as options.
%
% PAL_findMax:
% Modified: Palamedes version 1.2.0 (NP). Added 4D array. 
% Modified: Palamedes version 1.2.0 (NP). Modified such that routine works 
%   with array containing singleton dimensions. 
% Modified: Palamedes version 1.2.0 (NP). Fixed issue with function name 
%   (findMax -> PAL_findMax). 
% Modified: Palamedes version 1.2.0 (NP). Reduced memory load by
%   avoiding creation of maxVal array that existed in earlier version.
%
% PAL_Gumbel:
% Modified: Palamedes version 1.2.0 (NP). Added inverse PF and derivative 
%   of PF as options.
%
% PAL_HyperbolicSecant:
% Modified: Palamedes version 1.2.0 (NP). Added inverse PF and derivative 
%   of PF as options.
%
% PAL_inverseCumulativeNormal:
% Modified: Palamedes version 1.2.0 (NP). Added warning regarding removal 
%   of the function from a future version of Palamedes.
%
% PAL_inverseGumbel:
% Modified: Palamedes version 1.2.0 (NP). Added warning regarding removal 
%   of the function from a future version of Palamedes.
%
% PAL_inverseHyperbolicSecant:
% Modified: Palamedes version 1.2.0 (NP). Added warning regarding removal 
%   of the function from a future version of Palamedes.
%
% PAL_inverseLogistic:
% Modified: Palamedes version 1.2.0 (NP). Added warning regarding removal 
%   of the function from a future version of Palamedes.
%
% PAL_inverseWeibull:
% Modified: Palamedes version 1.2.0 (NP). Added warning regarding removal 
%   of the function from a future version of Palamedes.
%
% PAL_Logistic:
% Modified: Palamedes version 1.2.0 (NP). Added inverse PF and derivative 
%   of PF as options.
%
% PAL_MLDS_Bootstrap:
% Modified: Palamedes version 1.2.0 (NP): 'converged' is now array of 
%   logicals.
%
% PAL_PFLR_ModelComparison:
% Modified: Palamedes version 1.2.0 (NP): 'converged' is now array of 
%   logicals.
%
% PAL_PFML_BootstrapNonParametric:
% Modified: Palamedes version 1.2.0 (NP): 'converged' is now array of 
%   logicals.
% Modified: Palamedes version 1.2.0 (NP). Modified to accept 'searchGrid'
%   argument as a structure defining 4D parameter grid to search for
%   initial guesses for parameter values. See also
%   PAL_PFML_BruteForceFit.m.
%
% PAL_PFML_BootstrapNonParametricMultiple:
% Modified: Palamedes version 1.2.0 (NP): 'converged' is now array of 
%   logicals.
%
% PAL_PFML_BootstrapParametric:
% Modified: Palamedes version 1.2.0 (NP): 'converged' is now array of 
%   logicals.
% Modified: Palamedes version 1.2.0 (NP). Modified to accept 'searchGrid'
%   argument as a structure defining 4D parameter grid to search for
%   initial guesses for parameter values. See also
%   PAL_PFML_BruteForceFit.m.
%
% PAL_PFML_BootstrapParametricMultiple:
% Modified: Palamedes version 1.2.0 (NP): fixed error in help comments
%   (omission of two outputs in 'syntax' statement).
% Modified: Palamedes version 1.2.0 (NP): 'converged' is now array of 
%   logicals.
%
% PAL_PFML_Fit:
% Modified: Palamedes version 1.2.0 (NP). Modified to accept 'searchGrid'
%   argument as a structure defining 4D parameter grid to search for
%   initial guesses for parameter values. See also
%   PAL_PFML_BruteForceFit.m.
%
% PAL_PFML_GoodnessOfFit:
% Modified: Palamedes version 1.2.0 (NP): 'converged' is now array of 
%   logicals.
% Modified: Palamedes version 1.2.0 (NP). Modified to accept 'searchGrid'
%   argument as a structure defining 4D parameter grid to search for
%   initial guesses for parameter values. See also
%   PAL_PFML_BruteForceFit.m.
%
% PAL_PFML_GoodnessOfFitMultiple:
% Modified: Palamedes version 1.2.0 (NP): 'converged' is now array of 
%   logicals.
%
% PAL_PFML_GroupTrialsbyX:
% Modified: Palamedes version 1.2.0 (NP). Corrected error in function name
%   (Pal_GroupTrialsbyX -> PAL_PFML_GroupTrialsbyX)
%
% PAL_Scale0to1:
% Modified: Palamedes version 1.2.0 (NP): Modified to accept arrays of any
% size.
%
% PAL_Weibull:
% Modified: Palamedes version 1.2.0 (NP). Added inverse PF and derivative of
%   PF as options.
%
% Added routines:
%
% PAL_PFML_BruteForceFit:
% Introduced: Palamedes version 1.2.0 (NP): Fit PF using a brute-force 
%   search through 4D parameter space.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Version 1.3.0 Release September 19, 2011
%
% 
%
% Introduced PAL_PFML_lapseFit_Demo and PAL_PFML_gammaEQlambda_Demo to 
%   demonstrate new features.
%
% PAL_PF_SimulateObserverParametric;
% Modified: Palamedes version 1.3.0 (NP). Added options 'lapseFit' and
%   'gammaEQlambda'.
%
% PAL_PFLR_ModelComparison:
% Modified: Palamedes version 1.3.0 (NP). Added options 'lapseFit' and
%   'gammaEQlambda'.
%
% PAL_PFLR_TLR:
% Modified: Palamedes version 1.3.0 (NP). Added options 'lapseFit' and
%   'gammaEQlambda'.
%
% PAL_PFML_BootstrapNonParametric:
% Modified: Palamedes version 1.3.0 (NP). Added warning when 'LapseLimits'
%   argument is used but lapse is not a free parameter.
% Modified: Palamedes version 1.3.0 (NP). Added options 'lapseFit' and
%   'gammaEQlambda'.
%
% PAL_PFML_BootstrapNonParametricMultiple:
% Modified: Palamedes version 1.3.0 (NP). Added options 'lapseFit' and
%   'gammaEQlambda'.
%
% PAL_PFML_BootstrapParametric:
% Modified: Palamedes version 1.3.0 (NP). Added warning when 'LapseLimits'
%   argument is used but lapse is not a free parameter.
% Modified: Palamedes version 1.3.0 (NP). Added options 'lapseFit' and
%   'gammaEQlambda'.
%
% PAL_PFML_BootstrapParametricMultiple:
% Modified: Palamedes version 1.3.0 (NP). Added options 'lapseFit' and
%   'gammaEQlambda'.
%
% PAL_PFML_BruteForceFit:
% Modified: Palamedes version 1.3.0 (NP). Added options 'lapseFit' and
%   'gammaEQlambda'.
%
% PAL_PFML_DevianceGoF:
% Modified: Palamedes version 1.3.0 (NP). Added options 'lapseFit' and
%   'gammaEQlambda'.
%
% PAL_PFML_Fit:
% Modified: Palamedes version 1.3.0 (NP). Issue warning when 'LapseLimits'
%   argument is used but lapse is not a free parameter.
% Modified: Palamedes version 1.3.0 (NP). Added options 'lapseFit' and
%   'gammaEQlambda'.
%
% PAL_PFML_FitMultiple:
% Modified: Palamedes version 1.3.0 (NP). Added options 'lapseFit' and
%   'gammaEQlambda'.
%
% PAL_PFML_GoodnessOfFit:
% Modified: Palamedes version 1.3.0 (NP). Added warning when 'LapseLimits'
%   argument is used but lapse is not a free parameter.
% Modified: Palamedes version 1.3.0 (NP). Added options 'lapseFit' and
%   'gammaEQlambda'.
%
% PAL_PFML_GoodnessOfFitMultiple:
% Modified: Palamedes version 1.3.0 (NP). Added options 'lapseFit' and
%   'gammaEQlambda'.
%
% PAL_PFML_negLL:
% Modified: Palamedes version 1.3.0 (NP). Added options 'lapseFit' and
%   'gammaEQlambda'.
%
% PAL_PFML_negLLMultiple:
% Modified: Palamedes version 1.3.0 (NP). Fixed error in function name.
% Modified: Palamedes version 1.3.0 (NP). Added options 'lapseFit' and
%   'gammaEQlambda'.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Version 1.3.1 Release September 25, 2011
%
%Minor upgrade. Modified manner in which warnings regarding 'lapseFit' and
%   'lapseLimits' are issued. In Version 1.3.0 warnings would be issued by 
%   PAL_PFML_Fit if, for example, PAL_PFML_BootstrapParametric was called
%   without 'lapseFit' argument. All modifications related to above issue.
%
% PAL_PFML_negLL:
% Modified: Palamedes version 1.3.1 (NP). Added (hidden) option 'default' 
%   for 'lapseFit'.
%
% PAL_PFML_BruteForceFit:
% Modified: Palamedes version 1.3.1 (NP). Added (hidden) option 'default' 
%   for 'lapseFit'.
%
% PAL_PFML_Fit:
% Modified: Palamedes version 1.3.1 (NP). Added (hidden) option 'default' 
%   for 'lapseFit' and modified 'lapseFit' and 'lapseLimits' warnings to
%   avoid false throws of warnings.
%
% PAL_PFML_BootstrapParametric:
% Modified: Palamedes version 1.3.1 (NP). Added (hidden) option 'default' 
%   for 'lapseFit' and modified 'lapseFit' and 'lapseLimits' warnings to
%   avoid false throws of warnings.
%
% PAL_PFML_BootstrapNonParametric:
% Modified: Palamedes version 1.3.1 (NP). Added (hidden) option 'default' 
%   for 'lapseFit' and modified 'lapseFit' and 'lapseLimits' warnings to
%   avoid false throws of warnings.
%
% PAL_PFML_GoodnessOfFit:
% Modified: Palamedes version 1.3.1 (NP). Added (hidden) option 'default' 
%   for 'lapseFit' and modified 'lapseFit' and 'lapseLimits' warnings to
%   avoid false throws of warnings.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Version 1.4.0 Release February 15, 2012
%
% The major purpose of this upgrade is to make Palamedes compatible with
% GNU Octave. All Nelder-Mead simplex searches performed by the added
% routine PAL_minimize. Also added the option to constrain the guess rate
% in the PFML and PFLR functions. Some additional minor changes.
%
% PAL_AMPM_setupPM:
% Modified: Palamedes version 1.4.0 (NP): Fixed some stylistic nasties
%   ('Matlab style short circuit operator')
%
% PAL_AMRF_setupRF:
% Modified: Palamedes version 1.4.0 (NP): Fixed some stylistic nasties
%   ('Matlab style short circuit operator')
%
% PAL_AMRF_updateRF:
% Modified: Palamedes version 1.4.0 (NP): Fixed some stylistic nasties
%   ('Matlab style short circuit operator')
%
% PAL_AMUD_updateUD:
% Modified: Palamedes version 1.4.0 (NP): Fixed some stylistic nasties
%   ('Matlab style short circuit operator')
%
% PAL_CumulativeNormal:
% Modified: Palamedes version 1.4.0 (NP). Avoided use of erfcinv for
%   compatibility with Octave.
%
% PAL_inverseCumulativeNormal:
% Modified: Palamedes version 1.4.0 (NP). Functionality removed.
%
% PAL_inverseGumbel:
% Modified: Palamedes version 1.4.0 (NP). Functionality removed.
%
% PAL_inverseHyperbolicSecant:
% Modified: Palamedes version 1.4.0 (NP). Functionality removed.
%
% PAL_inverseLogistic:
% Modified: Palamedes version 1.4.0 (NP). Functionality removed.
%
% PAL_inverseWeibull:
% Modified: Palamedes version 1.4.0 (NP). Functionality removed.
%
% PAL_isIdentity:
% Modified: Palamedes version 1.4.0 (NP): Short-circuited logical 
%   operators.
%
% PAL_minimize:
% Introduced: Palamedes version 1.4.0 (NP)
%
% PAL_MLDS_Bootstrap:
% Modified: Palamedes version 1.4.0 (NP): Nelder-Mead simplex search
%   now performed by PAL_minimize.
%
% PAL_MLDS_Fit:
% Modified: Palamedes version 1.4.0 (NP): Nelder-Mead simplex search
%   now performed by PAL_minimize.
%
% PAL_PF_SimulateObserverParametric:
% Modified: Palamedes version 1.4.0 (NP). Short-circuited logical
%   operators.
%
% PAL_PFLR_ModelComparison:
% Modified: Palamedes version 1.4.0 (NP). Added 'guessLimits' option.
% Modified: Palamedes version 1.4.0 (NP). Short-circuited logical
%   operators.
%
% PAL_PFLR_TLR:
% Modified: Palamedes version 1.4.0 (NP). Added 'guessLimits' option.
% Modified: Palamedes version 1.4.0 (NP). Short-circuited logical 
%   operators.
%
% PAL_PFML_BootstrapNonParametric:
% Modified: Palamedes version 1.4.0 (NP). Added 'guessLimits' option.
%
% PAL_PFML_BootstrapNonParametricMultiple:
% Modified: Palamedes version 1.4.0 (NP). Check whether funcParamsSim needs
%   to be updated now performed properly.
% Modified: Palamedes version 1.4.0 (NP). Added 'guessLimits' option.
% Modified: Palamedes version 1.4.0 (NP). Short-circuited logical
%   operators.
%
% PAL_PFML_BootstrapParametric:
% Modified: Palamedes version 1.4.0 (NP). Added 'guessLimits' option.
% Modified: Palamedes version 1.4.0 (NP). Short-circuited logical
%   operators.
%
% PAL_PFML_BootstrapParametricMultiple:
% Modified: Palamedes version 1.4.0 (NP). Check whether funcParamsSim needs
%   to be updated now performed properly.
% Modified: Palamedes version 1.4.0 (NP). Added 'guessLimits' option.
% Modified: Palamedes version 1.4.0 (NP). Short-circuited logical
%   operators.
%
% PAL_PFML_Fit:
% Modified: Palamedes version 1.4.0 (NP). Added 'guessLimits' option.
% Modified: Palamedes version 1.4.0 (NP): Nelder-Mead simplex search
%   now performed by PAL_minimize.
% Modified: Palamedes version 1.4.0 (NP). Short-circuited logical
%   operators.
%
% PAL_PFML_FitMultiple:
% Modified: Palamedes version 1.4.0 (NP). Added 'guessLimits' option.
% Modified: Palamedes version 1.4.0 (NP): Nelder-Mead simplex search
%   now performed by PAL_minimize.
% Modified: Palamedes version 1.4.0 (NP). Short-circuited logical
%   operators.
%
% PAL_PFML_GoodnessOfFit:
% Modified: Palamedes version 1.4.0 (NP). Added 'guessLimits' option.
% Modified: Palamedes version 1.4.0 (NP). Short-circuited logical
%   operators.
%
% PAL_PFML_GoodnessOfFitMultiple:
% Modified: Palamedes version 1.4.0 (NP). Added 'guessLimits' option.
% Modified: Palamedes version 1.4.0 (NP). Short-circuited logical
%   operators.
%
% PAL_PFML_IndependentFit:
% Modified: Palamedes version 1.4.0 (NP): Short-circuited logical
%   operators.
%
% PAL_PFML_negLL:
% Modified: Palamedes version 1.4.0 (NP). Added 'guessLimits' option.
% Modified: Palamedes version 1.4.0 (NP). Short-circuited logical
%   operators.
%
% PAL_PFML_negLLMultiple:
% Modified: Palamedes version 1.4.0 (NP). Added 'guessLimits' option.
% Modified: Palamedes version 1.4.0 (NP). Short-circuited logical
%   operators.
%
% PAL_SDT_1AFCsameDiff_DiffMod_PHFtoDP:
% Modified: Palamedes version 1.4.0 (NP): Nelder-Mead simplex search now
%   performed by PAL-minimize.
%
% PAL_SDT_2AFCmatchSample_DiffMod_PCtoDP:
% Modified: Palamedes version 1.4.0 (NP): Nelder-Mead simplex search now
%   performed by PAL-minimize.
%
% PAL_SDT_2AFCmatchSample_DiffMod_PHFtoDP:
% Modified: Palamedes version 1.4.0 (NP): Nelder-Mead simplex search now
%   performed by PAL-minimize.
%
% PAL_SDT_2AFCmatchSample_IndMod_PCtoDP:
% Modified: Palamedes version 1.4.0 (NP): Nelder-Mead simplex search now
%   performed by PAL-minimize.
%
% PAL_SDT_2AFCmatchSample_IndMod_PHFtoDP:
% Modified: Palamedes version 1.4.0 (NP): Nelder-Mead simplex search now
%   performed by PAL-minimize.
%
% PAL_SDT_MAFC_DPtoPC:
% Modified: Palamedes version 1.4.0 (NP): Solved Octave incompatibility
%   issue.
%
% PAL_SDT_MAFC_PCtoDP:
% Modified: Palamedes version 1.4.0 (NP): Nelder-Mead simplex search now
%   performed by PAL-minimize.
%
% PAL_SDT_MAFCmatchSample_DiffMod_PCtoDP:
% Modified: Palamedes version 1.4.0 (NP): Nelder-Mead simplex search now
%   performed by PAL-minimize.
%
% PAL_SDT_MAFCoddity_PCtoDP:
% Modified: Palamedes version 1.4.0 (NP): Nelder-Mead simplex search now
%   performed by PAL_minimize.
%
% PAL_spreadPF:
% Modified: Palamedes version 1.4.0 (NP) Avoided use of erfcinv for
%   compatibility with Octave.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Version 1.4.1 Release February 18, 2012
%
% Fixed bug which resulted when PAL_minimize encounters NaN's: 
%   PAL_PFML_negLL no longer returns NaN's.
%
% PAL_PFML_negLL:
% Modified: Palamedes version 1.4.1 (NP). Ignore NaN's (which arise when
%   0*log(0) is evaluated).
% 
% PAL_nansum:
% Introduced: Palamedes version 1.4.1 (NP).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Version 1.4.2 Release March 17, 2012
%
% Fixes a bug that would make PAL_PFML_BootstrapParametricMultiple and
%   PAL_PFML_FitMultiple misbehave when data arrays (StimLevels,
%   NumPos, OutOfNum) were such that a call to PAL_PFML_GroupTrialsbyX 
%   would have changed their arrangement.
%
% PAL_PFML_BootstrapParametricMultiple:
% Modified: Palamedes version 1.4.2 (NP). Fixed typo-bug Stimlevels -
%   StimLevels (line 207 in version 1.4.2). This bug would have caused
%   trouble if this function was used with data-arrays that would have 
%   changed by running it through PAL_PFML_GroupTrialsbyX
%
% PAL_PFML_FitMultiple:
% Modified: Palamedes version 1.4.2 (NP). Fixed typo-bug Stimlevels ->
%   StimLevels (line 330 in version 1.4.2). This bug would have caused
%   trouble if this function was used with data-arrays that would have 
%   changed by running it through PAL_PFML_GroupTrialsbyX
%
% PAL_PFML_negLL: 
% Modified: Palamedes version 1.4.2 (NP). Changed manner in which NaN's
%   arising from evaluating 0*log(0) are handled (some more).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Version 1.4.3 Release March 19, 2012
%
% Undoes an unintended change introduced in 1.4.2. In PAL_minimize, the
%   value of delta was changed to 0.02. Now it's back to 0.05.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Version 1.4.4 Release March 29, 2012
%
% Undoes another unintended change introduced in 1.4.2. In PAL_PFLR_TLR, 
%   procedure was 'paused' in case convergence failure occurred (for 
%   debugging purposes). Now undone.
%
%PAL_CumulativeNormal: 
%Modified: Palamedes version 1.4.4 (NP). Fixed bug: '.5*erfc' -> '.5.*erfc'
%   credit: Loes van Dam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Version 1.4.5 Release April 6, 2012
%
% PAL_AMUD_setupUD:
% Modified: Palamedes version 1.4.5 (NP) Initialized UD.xMax and UD.xMin to
%   Inf and -Inf, respectively.
% Modified: Palamedes version 1.4.5 (NP): Bug fix: Changed 'UD.truncate ==
%   1' to 'strcmp(UD.truncate,'yes')' (line 68 in 1.4.5)
% Modified: Palamedes version 1.4.5 (NP): Bug fix: Changed 'UD.reversal < 1'
%   to max(UD.reversal) < 1 (lines 47 & 66 in 1.4.5)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Version 1.5.0 Release June 10, 2012
%
% Purpose of release: Implementation of Psi+ method. The Psi+ method
% modifies the Psi method to allow optimization of guess rate and/or lapse
% rate as well as threshold and slope.
%
% A minor incompatibility with older versions is introduced. See warning in
% PAL_AMPM_setupPM or www.palamedestoolbox.org/pal_ampm_incompatibility.html
%
% PAL_AMPM_createLUT:
% Modified: Palamedes version 1.5.0 (NP): Modified to allow 'Psi+' 
%   (optimization of guess rate and/or lapse rate as well as threshold and 
%   slope).
%
% PAL_AMPM_posteriorTplus1:
% Modified: Palamedes version 1.5.0 (NP): Modified to allow 'Psi+' 
%   (optimization of guess rate and/or lapse rate as well as threshold and 
%   slope).
%
% PAL_AMPM_pSuccessGivenx:
% Modified: Palamedes version 1.5.0 (NP): Modified to allow 'Psi+' 
%   (optimization of guess rate and/or lapse rate as well as threshold and 
%   slope).
%
% PAL_AMPM_setupPM:
% Modified: Palamedes version 1.5.0 (NP): Modified to allow 'Psi+' 
%   (optimization of guess rate and/or lapse rate as well as threshold and 
%   slope).
%
% PAL_AMPM_updatePM:
% Modified: Palamedes version 1.5.0 (NP): Modified to allow 'Psi+' 
%   (optimization of guess rate and/or lapse rate as well as threshold and 
%   slope).
%
% PAL_Entropy:
% Modified: Palamedes version 1.5.0 (NP): Modified to accept N-D arrays.
%   Added option to return entropies across limited number of dimensions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Version 1.6.0 Release March 15, 2013
%
% Purposes of release: (1) Introduce new SDT routines: routines that fit 
%   ROC curve, determine standard errors of its parameters, goodness of the 
%   fit, and perform a model comparison to determine whether the ratio of 
%   standard deviations of signal and noise distributions equals 1 (or any 
%   other value). Also introduced 3AFC and MAFC oddity routines (2) Extend 
%   functionality of Psi method routines: when lapse rate is included in 
%   the posterior distribution it can now be assumed that gamma and lambda 
%   are equal (e.g., in bistable percept task in which both gamma and 
%   lambda can be thought of as estimating the lapse rate). User can also 
%   marginalize nuisance parameters such that expected entropy is 
%   determined across parameters of interest only.
%
% PAL_SDT_cumulateHF:
% Introduced: Palamedes version 1.6.0 (FK & NP)
%
% PAL_SDT_3AFCoddity_IndMod_DPtoPC
% Introduced: Palamedes version 1.6.0 (FK)
%
% PAL_SDT_3AFCoddity_IndMod_DPtoPCpartFuncA
% Introduced: Palamedes version 1.6.0 (FK)
%
% PAL_SDT_3AFCoddity_IndMod_DPtoPCpartFuncB
% Introduced: Palamedes version 1.6.0 (FK)
%
% PAL_SDT_3AFCoddity_IndMod_PCtoDP
% Introduced: Palamedes version 1.6.0 (FK)
%
% PAL_SDT_MAFCoddity_IndMod_DPtoPC
% Introduced: Palamedes version 1.6.0 (FK)
%
% PAL_SDT_MAFCoddity_IndMod_PCtoDP
% Introduced: Palamedes version 1.6.0 (FK)
%
% PAL_SDT_ROC_SimulateObserverParametric:
% Introduced: Palamedes version 1.6.0 (FK & NP)
%
% PAL_SDT_ROCML_BootstrapParametric:
% Introduced: Palamedes version 1.6.0 (FK & NP)
%
% PAL_SDT_ROCML_Fit:
% Introduced: Palamedes version 1.6.0 (FK & NP)
%
% PAL_SDT_ROCML_GoodnessOfFit:
% Introduced: Palamedes version 1.6.0 (FK & NP)
%
% PAL_SDT_ROCML_negLL:
% Introduced: Palamedes version 1.6.0 (FK & NP)
%
% PAL_SDT_ROCML_negLLNonParametric:
% Introduced: Palamedes version 1.6.0 (FK & NP)
%
% PAL_SDT_ROCML_RatioSDcomparison:
% Introduced: Palamedes version 1.6.0 (FK & NP)
%
% PAL_Quick:
% Introduced: Palamedes version 1.6.0 (NP)
%
% PAL_logQuick:
% Introduced: Palamedes version 1.6.0 (NP)
%
% PAL_AMPM_CreateLUT:
% Modified: Palamedes version 1.6.0 (NP): Modified to allow Psi method to
% 	assume that gamma equals lambda.
%
% PAL_AMPM_setupPM:
% Modified: Palamedes version 1.6.0 (NP): Modified to allow Psi method to
%   assume that gamma equals lambda (e.g., both estimating lapse rate in a
%   bistable percept task). Modified to allow marginalization of
%   parameters before entropy is calculated.
%
% PAL_AMPM_updatePM:
% Modified: Palamedes version 1.6.0 (NP): Modified to allow Psi method to
%   assume that gamma equals lambda (e.g., both estimating lapse rate in a
%   bistable percept task). Modified to allow marginalization of
%   parameters before entropy is calculated.
%
% PAL_Entropy:
% Modified: Palamedes version 1.6.0 (NP): Modified to allow marginalizing
%   out dimension(s) before calculating entropy.
%
% PAL_findMax:
% Modified: Palamedes version 1.6.0 (NP): Revamped (and simplified) 
%   strategy in order to allow finding maximum in N-D arrays.
%
% PAL_PFML_bruteForceFit:
% Modified: Palamedes version 1.6.0 (NP): Avoid waste of memory (in case
%   user sets 'gammaEQlambda' to 'true' but supplies non-scalar
%   searchGrid.gamma) by setting searchGrid.gamma to scalar with 
%   (arbitrary) value of 0.
%
% PAL_PFML_negLLNonParametric:
% Modified: Palamedes version 1.6.0 (NP): Deal with 0.*log(0) issue using
%   PAL_nansum.
%
% PAL_PFML_TtoP:
% Modified: Palamedes version 1.6.0 (NP): Allow one to fix a parameter in
%   one condition while allowing it to be free in other conditions when 
%   using a contrast matrix to constrain parameter in functions that fit 
%   PFs simultaneously to multiple conditions. Example: passing 
%   [1 0 0; 0 1 0] to, say, 'thresholds' in PAL_PFML_FitMultiple fixes the
%   threshold value in the third condition while estimating independent 
%   thresholds in conditions 1 and 2. Example 2: [0 1 1 1; 0 -1 0 1] will
%   fix parameter in condition 1 while fitting an intercept and linear 
%   trend to conditions 2 through 4.
%
% PAL_SDT_MAFC_DPtoPCpartFunc:
% Modified: Palamedes version 1.6.0 (NP): Modified to eliminate reliance on
%   Matlab Statistics Toolbox.
%
% PAL_SDT_1AFC_DPtoPHF:
% Modified: Palamedes version 1.6.0 (FK & NP): Modified to allow
%   functionality of new ROC routines.
%
% PAL_SDT_1AFC_PHFtoDP:
% Modified: Palamedes version 1.6.0 (FK & NP): Modified to allow
%   functionality of new ROC routines.
%
% PAL_PFLR_ModelComparison:
% Modified: Palamedes version 1.6.0 (NP): Use of 'lapseLimits' option now
%   allowed when lapse rates in model(s) defined using model matrix or
%   custom reparametrization.
%
% PAL_PFML_BootstrapNonParametricMultiple:
% Modified: Palamedes version 1.6.0 (NP): Use of 'lapseLimits' option now
%   allowed when lapse rates in model(s) defined using model matrix or
%   custom reparametrization.
%
% PAL_PFML_BootstrapParametricMultiple:
% Modified: Palamedes version 1.6.0 (NP): Use of 'lapseLimits' option now
%   allowed when lapse rates in model(s) defined using model matrix or
%   custom reparametrization.
%
% PAL_PFML_FitMultiple:
% Modified: Palamedes version 1.6.0 (NP): Use of 'lapseLimits' option now
%   allowed when lapse rates in model(s) defined using model matrix or
%   custom reparametrization.
%
% PAL_PFML_GoodnessOfFitMultiple:
% Modified: Palamedes version 1.6.0 (NP): Use of 'lapseLimits' option now
%   allowed when lapse rates in model(s) defined using model matrix or
%   custom reparametrization.
%
% PAL_PFML_negLLMultiple:
% Modified: Palamedes version 1.6.0 (NP): Modified to allow modification
%   made to PAL_PFLR_ModelComparison, 
%   PAL_PFML_BootstrapNonParametricMultiple,
%   PAL_PFML_BootstrapParametricMultiple, PAL_PFML_FitMultiple, and
%   PAL_PFML_GoodnessOfFitMultiple.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Version 1.6.1 Release December 20, 2013
%
% Purposes of release: (1) Reduction of computational load of AMPM 
% routines: duplication of a computation avoided by absorbing 
% PAL_AMPM_pSuccessGivenx function into PAL_AMPM_PosteriorTplus1, 
% PAL_AMPM_PosteriorTplus1 recoded to increase efficiency. 
% Functionality of user-end AMPM routines not affected (2) update 
% PAL_spreadPF (3) fix incompatibility of PAL_AMPM_DEMO with
% GNU Octave.
%
% PAL_AMPM_pSuccessGivenx:
% Removed: Palamedes version 1.6.1 (NP): Removed. Functionality absorbed in
%   PAL_AMPM_PosteriorTplus1 in order to reduce computational load of AMPM 
%   routines.
%
% PAL_AMPM_expectedEntropy:
% Modified: Palamedes version 1.6.1 (NP): pSuccessGivenx imported from
%   PAL_AMPM_PosteriorTSplus1 rather than from (the now removed) 
%   PAL_AMPM_pSuccessGivenx.
%
% PAL_AMPM_PosteriorTplus1:
% Modified: Palamedes version 1.6.1 (NP): Recoded to increase efficiency, 
%   absorbed functionality of (the now removed) PAL_AMPM_pSuccessGivenx.
%
% PAL_AMPM_setupPM:
% Modified: Palamedes version 1.6.1 (NP): Modified to comply with
%   absorption of PAL_AMPM_pSuccessGivenx functionality into
%   PAL_AMPM_PosteriorTplus1.
%
% PAL_spreadPF:
% Modified: Palamedes version 1.6.1 (NP): Added the options to compute
%   spread for Quick and logQuick functions (both introduced in version
%   1.6.0).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Version 1.6.2 Release February 17, 2014
%
% Purposes of release: Improve PAL_Contrast (functionality and
%   user-friendliness), minor maintenance on PAL_minimize.
% 
% PAL_isOrthogonal:
% Introduced: Palamedes version 1.6.2 (NP): Checks orthogonality of the
%   pairings of row vectors. Called by PAL_Contrasts to verify 
%   orthogonality of model matrices.
%
% PAL_Contrasts:
% Modified: Palamedes version 1.6.2 (NP): Modified to return integer-
%   valued contrast coefficients when possible, also verifies orthogonality
%   of polynomial contrasts and issues warning when incapable of finding a
%   full set of orthogonal polynomial contrasts, allows user to specify IV
%   levels that are not equally spaced.
%
% PAL_minimize:
% Modified: Palamedes version 1.6.2 (NP): Directs user to possible cause of
%   issued warning regarding values of MaxFunEvals and/or MaxIter.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Version 1.6.3 Release March 3, 2014
%
% Purposes of release: Restore backward compatibility to Matlab Release 14
%   SP1 (compatibility checked with Matlab version 7.0.1.24704, bare Matlab
%   [no Matlab toolboxes]) and restore compatibility with GNU Octave 
%   (compatibility checked with GNU Octave version 3.8.0 under Ubuntu 
%   13.10). Minor maintenance.
%
% PAL_AMPM_PosteriorTplus1:
% Modified: Palamedes version 1.6.3 (NP): Checks for presence of bsxfun.m
%   and performs alternative operation if not present (in order to restore 
%   compatibility with older versions of Matlab).
%
% PAL_Contrasts:
% Modified: Palamedes version 1.6.3 (NP): No longer changes type of 'temp'
%   to int64. (in order to restore compatibility with older versions of 
%   Matlab).
%
% PAL_PFML_BootstrapNonParametric:
% Modified: Palamedes version 1.6.3 (NP): No longer calls iscolumn.m (in 
%   order to restore compatibility with older versions of Matlab).
%
% PAL_PFML_BootstrapParametric:
% Modified: Palamedes version 1.6.3 (NP): No longer calls iscolumn.m (in 
%   order to restore compatibility with older versions of Matlab).
%
% PAL_PFML_GoodnessOfFit:
% Modified: Palamedes version 1.6.3 (NP): No longer calls iscolumn.m (in 
%   order to restore compatibility with older versions of Matlab).
%
% PAL_PFML_PtoT:
% Modified: Palamedes version 1.6.3 (NP): Removed rounding of theta values 
%   that are very near zero. Problem is taken care of in PAL_minimize.
%
% PAL_Quick:
% Modified: Palamedes version 1.6.3 (NP): Negative stimulus intensity
%   values will generate an error. Error message suggests to use
%   PAL_logQuick.
%
% PAL_Weibull:
% Modified: Palamedes version 1.6.3 (NP): Negative stimulus intensity
%   values will generate an error. Error message suggests to use
%   PAL_logQuick.
% 
% PAL_AMPM_CreateLUT and
% PAL_AMPM_expectedEntropy and
% PAL_AMPM_PosteriorTplus1 and
% PAL_AMPM_setupPM and
% PAL_AMPM_updatePM and
% PAL_AMRF_pdfDescriptives and
% PAL_AMRF_setupRF and
% PAL_AMRF_updateRF and
% PAL_AMUD_analyzeUD and
% PAL_AMUD_setupUD and
% PAL_Contrasts and
% PAL_CumulativeNormal and
% PAL_Entropy and
% PAL_findMax and
% PAL_Gumbel and
% PAL_HyperbolicSecant and
% PAL_isOrthogonal and
% PAL_Logistic and
% PAL_logQuick and
% PAL_MeanSDSSandSE and
% PAL_minimize and
% PAL_MLDS_Bootstrap and
% PAL_MLDS_Fit and
% PAL_MLDS_GroupTrialsbyX and
% PAL_MLDS_SimulateObserver and
% PAL_PF_SimulateObserverParametric and
% PAL_PFBA_Fit and
% PAL_PFLR_ModelComparison and
% PAL_PFLR_setupMC and
% PAL_PFLR_TLR and
% PAL_PFML_BootstrapNonParametric and
% PAL_PFML_BootstrapNonParametricMultiple and
% PAL_PFML_BootstrapParametric and
% PAL_PFML_BootstrapParametricMultiple and
% PAL_PFML_BruteForceFit and
% PAL_PFML_DevianceGoF and
% PAL_PFML_Fit and
% PAL_PFML_FitMultiple and
% PAL_PFML_GoodnessOfFit and
% PAL_PFML_GoodnessOfFitMultiple and
% PAL_PFML_GroupTrialsbyX and
% PAL_PFML_LLsaturated and
% PAL_PFML_negLL and
% PAL_PFML_negLLMultiple and
% PAL_PFML_negLLNonParametric and
% PAL_PFML_paramsTry and
% PAL_PFML_PtoT and
% PAL_Quick and
% PAL_randomizeArray and
% PAL_SDT_1AFC_DPtoPHF and
% PAL_SDT_1AFC_PHFtoDP and
% PAL_SDT_1AFCsameDiff_DiffMod_PHFtoDP and
% PAL_SDT_1AFCsameDiff_IndMod_PHFtoDP and
% PAL_SDT_2AFC_PHFtoDP and
% PAL_SDT_2AFCmatchSample_DiffMod_PCtoDP and
% PAL_SDT_2AFCmatchSample_DiffMod_PHFtoDP and
% PAL_SDT_2AFCmatchSample_IndMod_PCtoDP and
% PAL_SDT_2AFCmatchSample_IndMod_PHFtoDP and
% PAL_SDT_3AFCoddity_IndMod_DPtoPC and
% PAL_SDT_3AFCoddity_IndMod_PCtoDP and
% PAL_SDT_MAFC_DPtoPC and
% PAL_SDT_MAFC_PCtoDP and
% PAL_SDT_MAFCmatchSample_DiffMod_DPtoPC and
% PAL_SDT_MAFCmatchSample_DiffMod_PCtoDP and
% PAL_SDT_MAFCoddity_DPtoPC and
% PAL_SDT_MAFCoddity_IndMod_DPtoPC and
% PAL_SDT_MAFCoddity_IndMod_PCtoDP and
% PAL_SDT_MAFCoddity_PCtoDP and
% PAL_SDT_ROCML_BootstrapParametric and
% PAL_SDT_ROCML_Fit and
% PAL_SDT_ROCML_GoodnessOfFit and
% PAL_SDT_ROCML_negLL and
% PAL_SDT_ROCML_negLLNonParametric and
% PAL_SDT_ROCML_RatioSDcomparison and
% PAL_spreadPF and
% PAL_unpackParamsPF and
% PAL_Weibull:
% Modified: Palamedes version 1.6.3 (NP): Minor maintenance: Removed
%   suppression of logOfZero warnings, removed suppression of DivideByZero 
%   warnings (users of older versions of Matlab and users of Octave may
%   occassionally encounter either warning, please ignore [they are 
%   inconsequential] or suppress in code), Palamedes warnings are given
%   identifiers, separated output variables using commas, pre-allocated
%   variables, etc.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version 1.6.3.1 Release April 10, 2014
%
% Purposes of release: Direct user to page:
% www.palamedestoolbox.org/weibullandfriends.html in help comments of:
% PAL_Gumbel, PAL_logQuick, PAL_Quick, and PAL_Weibull. Error messages
% issued by PAL_Weibull and PAL_Quick when negative stimulus values are
% passed direct user to same page.
%
% No change in functionality.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version 1.7.0 Release June 27, 2014
%
% Purposes of release: (1) Introduce SDT summation routines (listed below)
% (2) Introduce PAL_LLfixed and PAL_GoodnessOfFitZeroDF (3) Correct some 
% regressions that occurred in 1.6.3.1 (some of the updates introduced in 
% version 1.6.3 were, somehow, undone in 1.6.3.1). Version 1.7.0 was 
% created by starting from 1.6.3, performing upgrade to 1.6.3.1 again, then 
% making the changes that are recorded below to result in version 1.7.0.
%
% PAL_SDT_AS_PCtoSL:
% Introduced: Palamedes version 1.7.0 (FK&NP)
% uses an iterative search procedure of PAL_SDT_AS_SLtoPC to calculate the 
% stimulus level x from proportion correct PC, a PC that would result 
% from the ADDITIVE SUMMATION of n identical stimuli in Q monitored 
% channels, if x was subject to a scaling factor g and a transducer 
% exponent p such that d'= (gx)^p, for a M-AFC 
% (M-alternative-forced-choice) task under the assumptions of 
% signal-detection-theory, assuming an unbiased observer.
%
% PAL_SDT_AS_uneqSLtoPC:
% Introduced: Palamedes version 1.7.0 (FK&NP)
% uses numerical integration to calculate proportion correct PC based 
% on the ADDITIVE SUMMATION of unequal stimulus levels x in Q monitored 
% channels, with each stimulus level x subject to a scaling factor g 
% and a transducer exponent p such that d'= (gx)^p, for a M-AFC 
% (M-alternative-forced-choice) task under the assumptions of SDT and 
% assuming an unbiased observer
%
% PAL_SDT_PS_DPtoPCpartFunc:
% Introduced: Palamedes version 1.7.0 (FK&NP)
% Internal function
%
% PAL_SDT_PS_MonteCarlo_SLtoPC:
% Introduced: Palamedes version 1.7.0 (FK&NP)
% uses numerical integration to calculate proportion correct PC based 
% on the PROBABILITY SUMMATION of n identical stimuli in Q monitored 
% channels, with signal level x subject to a scaling factor g and 
% transducer exponent p such that d'= (gx)^p, for a M-AFC 
% (M-alternative-forced-choice) task under the assumptions of 
% signal-detection-theory and assuming an unbiased observer
%
% PAL_SDT_PS_MonteCarlo_uneqSLtoPC:
% Introduced: Palamedes version 1.7.0 (FK&NP)
% uses Monte Carlo simulation to calculate proportion correct based 
% on the PROBABILITY SUMMATION of unequal stimulus levels x in Q monitored 
% channels, with each stimulus level x subject to a scaling factor g 
% and a transducer exponent p such that d'= (gx)^p, for a M-AFC 
% (M-alternative-forced-choice) task, under the assumptions of SDT and 
% assuming an unbiased observer
%
% PAL_SDT_PS_PCtoSL:
% Introduced: Palamedes version 1.7.0 (FK&NP)
% uses an iterative search procedure of PAL_SDT_AS_SLtoPC to calculate 
% stimulus level x from proportion correct PC, a PC that would result 
% from the PROBABILITY SUMMATION of n identical stimuli in Q monitored 
% channels, with x subject to a scaling factor g and transducer exponent p 
% such that d'= (gx)^p, for a M-AFC (M-alternative-forced-choice) task 
% under the assumptions of signal-detection-theory, assuming an unbiased 
% observer.
%
% PAL_SDT_PS_SLtoPC:
% Introduced: Palamedes version 1.7.0 (FK&NP)
% uses numerical integration to calculate proportion correct PC based 
% on the PROBABILITY SUMMATION of n identical stimuli in Q monitored 
% channels, with stimulus level x subject to a scaling factor g and 
% transducer exponent p such that d'= (gx)^p, for a M-AFC 
% (M-alternative-forced-choice) task, under the assumptions of 
% signal-detection-theory and assuming an unbiased observer
%
% PAL_SDT_PS_uneqDPtoPCpartFunc:
% Introduced: Palamedes version 1.7.0 (FK&NP)
% Internal function
%
% PAL_SDT_PS_uneqDPtoPCpartFunc2:
% Introduced: Palamedes version 1.7.0 (FK&NP)
% Internal function
%
% PAL_SDT_PS_uneqSLtoPC:
% Introduced: Palamedes version 1.7.0 (FK&NP)
% uses numerical integration to calculate proportion correct based 
% on the PROBABILITY SUMMATION of unequal stimulus levels in Q monitored 
% channels, with each stimulus level x subject to a scaling factor g 
% and a transducer with exponent p such that d'= (gx)^p, for a M-AFC 
% (M-alternative-forced-choice) task, under the assumptions of SDT and 
% assuming an unbiased observer
%
% PAL_LLfixed:
% Introduced: Palamedes version 1.7.0 (NP)
% Internal Function
%
% PAL_GoodnessOfFitZeroDF:
% Introduced: Palamedes version 1.7.0 (NP)
% Determine Goodness-of-Fit of model of binomial data where the model 
% consists only of fixed parameter values (i.e., the model has zero degrees 
% of freedom).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version 1.8.0 Release February 19, 2015
%
% Purpose of release is to a) modify the names of the SDT (signal-
% detection-theory) differencing model oddity routines, b) introduce 2AFC 
% SDT routines, c) introduce routines for fitting PFs (psychometric 
% functions) with SDT models, d) introduce new routines for summation under 
% SDT and e) introduce routines for fitting summation PFs under SDT.
% 
% PAL_SDT_MAFC_Oddity_DPtoPC
% Modified: Palamedes version 1.8.0 (FK)
% Disabled, renamed PAL_SDT_MAFC_Oddity_DiffMod_DPtoPC
% 
% PAL_SDT_MAFC_Oddity_DPtoPC
% Modified: Palamedes version 1.8.0 (FK)
% Disabled, renamed PAL_SDT_MAFC_Oddity_DiffMod_PCtoDP
% 
% PAL_SDT_MAFC_Oddity_DiffMod_DPtoPC
% Introduced: Palamedes version 1.8.0 (FK & NP)
% PAL_SDT_MAFCoddity_DiffMod_DPtoPC converts d'(d-prime) to 
% proportion correct for a M-AFC (M-alternative-forced-choice) 
% oddity task, assuming a Differencing model and an unbiased 
% observer
% 
% PAL_SDT_MAFC_Oddity_DiffMod_PCtoDP
% Introduced: Palamedes version 1.8.0 (FK & NP)
% PAL_SDT_MAFCoddity_DiffMod_PCtoDP converts proportion correct into 
% d'(d-prime) for a M-AFC (M-alternative-forced-choice) oddity
% task, assuming a Differencing Observer model and an unbiased observer
% 
% PAL_SDT_2AFC_DPtoPC
% Introduced: Palamedes version 1.8.0 (FK)
% PAL_SDT_2AFC_DPtoPC converts d'(d-prime) into proportion correct for a 
% 2AFC (two-alternative-forced-choice) task
% 
% PAL_SDT_2AFC_PCtoDP
% Introduced: Palamedes version 1.8.0 (FK)
% PAL_SDT_2AFC_PCtoDP converts proportion correct to d' (d-prime) for a 
% 2AFC (two-alternative-forced-choice) task
%
% PAL_SDT_PS_2uneqSLtoPC
% Introduced: Palamedes version 1.8.0 (FK)
% uses numerical integration to calculate proportion correct based 
% on the PROBABILITY SUMMATION of two unequal intensity stimuli, one 
% with an intensity of x the other r*x, in Q monitored channels, with each 
% x subject to a scaling factor g and transducer exponent p 
% such that d'= (gx)^p, in a M-AFC (M-alternative-forced-choice) task, 
% under the assumptions of SDT and assuming an unbiased observer
% 
% PAL_SDT_AS_2uneqSLtoPC
% Introduced: Palamedes version 1.8.0 (FK)
% uses numerical integration to calculate proportion correct based 
% on the ADDITIVE SUMMATION of two unequal stimulus intensities, one of 
% which is x the other r*x, in Q monitored channels, with each x subject 
% to a scaling factor g and transducer exponent p such that d'= (gx)^p, 
% for a M-AFC (M-alternative-forced-choice) task, under the assumptions of 
% SDT and assuming an unbiased observer
% 
% PAL_SDT_PS_PCto2uneqSL
% Introduced: Palamedes version 1.8.0 (FK)
% uses an iterative search procedure of PAL_SDT_PS_2uneqSLtoPC to 
% calculate two stimulus levels x with a ratio r from proportion correct PC, 
% a PC that would result from the PROBABILITY SUMMATION of the two 
% different stimuli, in Q monitored channels, with x subject to a scaling 
% factor g and transducer exponent p such that d'= (gx)^p, for a M-AFC 
% (M-alternative-forced-choice) task under the assumptions of 
% signal-detection-theory, assuming an unbiased observer.
% 
% PAL_SDT_AS_PCto2uneqSL
% Introduced: Palamedes version 1.8.0 (FK)
% uses an iterative search procedure of PAL_SDT_AS_2uneqSLtoPC to 
% calculate two stimulus levels x with ratio r from proportion correct PC, 
% a PC that would result from the ADDITIVE SUMMATION of the two 
% x values, in Q monitored channels, with x subject to a scaling 
% factor g and transducer exponent p such that d'= (gx)^p, for a M-AFC 
% (M-alternative-forced-choice) task under the assumptions of 
% signal-detection-theory, assuming an unbiased observer.
% 
% PAL_SDT_SLtoPC
% Introduced: Palamedes version 1.8.0 (FK)
% Converts stimulus intensities to proportion correct assuming parameters 
% of an SDT (signal detection theory) model
% 
% PAL_SDT_PCtoSL
% Introduced: Palamedes version 1.8.0 (FK)
% Converts proportion correct to stimulus intensity assuming parameters of 
% an SDT (signal detection theory) model
% 
% PAL_SDT_PFML_negLL
% Introduced: Palamedes version 1.8.0 (FK & NP)
% (negative) Log Likelihood associated with fit of SDT (signal-detection-
% theory) psychometric function 
% 
% PAL_SDT_PFML_Fit
% Introduced: Palamedes version 1.8.0 (FK & NP)
% Fits the parameters of an SDT (signal dtection theory) psychometric 
% function (PF) using a Maximum Likelihood criterion.
% 
% PAL_SDT_PF_SimulateObserverParametric
% Introduced: Palamedes version 1.8.0 (FK & NP)
% Simulate observer characterized by a PF (psychomtric function) from an 
% SDT (signal detection thory) model.
% 
% PAL_SDT_PFML_BootstrapParametric
% Introduced: Palamedes version 1.8.0 (FK & NP)
% Perform parametric bootstrap to determine standard errors on parameters 
% of an SDT (signal detection theory) model fitted to PF (psychometric 
% function) data
% 
% PAL_SDT_PFML_GoodnessOfFit
% Introduced: Palamedes version 1.8.0 (FK & NP)
% Determines goodness-of-fit of the parameters of a signal detection theory 
% (SDT) model fitted to psychometric function (PF) data
% 
% PAL_SDT_Summ_MultiplePFML_negLL
% Introduced: Palamedes version 1.8.0 (FK & NP)
% (negative) Log Likelihood associated with fit of summation psychometric 
% function 
% 
% PAL_SDT_Summ_MultiplePFML_Fit
% Introduced: Palamedes version 1.8.0 (FK & NP)
% Fits the parameters of a set of summation psychometric curves using a 
% Maximum Likelihood criterion, based on either a probability or additive 
% summation model
% 
% PAL_SDT_Summ_MultiplePF_SimulateObserverParametric
% Introduced: Palamedes version 1.8.0 (FK & NP)
% Simulate observer for a summation experiment with multiple PFs 
% (psychometric functions) based on SDT (signal detection theory)
% 
% PAL_SDT_Summ_MultiplePFML_BootstrapParametric
% Introduced: Palamedes version 1.8.0 (FK & NP)
% Perform parametric bootstrap to determine standard errors on parameters 
% of SDT (signal detection theory) probability and additive summation 
% models simultaneously fitted to multiple PF (psychometric function) data
% 
% PAL_SDT_Summ_MultiplePFML_GoodnessOfFit
% Introduced: Palamedes version 1.8.0 (FK & NP)
% Obtains Goodness-of-fit of parameters of SDT (signal detection theory) 
% probability and additive summation models obtained from simultaneously 
% fitted multiple PF (psychometric function) data  
