%   
%   PAL_PFML_CustomDefine  Instructions on the use of custom
%   reparametrization of parameters in PF fitting.
%
%   This information applies to custom reparametrization of parameters in 
%   the functions:
%       PAL_PFML_FitMultiple
%       PAL_PFML_BootstrapParametricMultiple
%       PAL_PFML_BootstrapNonParametricMultiple
%       PAL_PFML_GoodnessOfFitMultiple
%       PAL_PFLR_ModelComparison
%
%   In the help comments of the above functions it is explained how the 
%   user can constrain the parameters of psychometric functions between
%   several data sets using the 'fixed', 'constrained', or 'unconstrained'
%   arguments. For example, the call:
%
%   paramsFitted = PAL_PFML_FitMultiple(StimLevels,NumPos,OutOfNum,...
%       params,PF,'thresholds','unconstrained','slopes','constrained',...
%       'guessrates','fixed','lapserates','fixed');
%
%   instructs Palamedes to fit each data set with individual thresholds,
%   but to fit a single shared value for the slope, and to use fixed
%   values for the guessrate and lapse rate. Also explained there is how
%   the user may reparametrize any of the PFs parameters using linear 
%   contrasts. Here it is explained how the user may custom-define non-
%   linear constraints on the PF parameters. It is assumed here that the 
%   user is already familiar with the general usage of PAL_PFML_FitMultiple. 
%   Type 'help PAL_PFML_FitMultiple' for more information.
%
%   We explain by way of example.
%
%   Imagine one is interested in modeling the improvement in performance on 
%   some task with practice. Such improvement might be exhibited by a 
%   decreasing threshold across consecutive sessions. Such improvements may 
%   be modeled/summarized by a learning curve, a common parametrization of 
%   which is the 'exponential decay function':
%
%   Threshold(session) = a + b*exp(-c*(session-1))
%
%   Parameter 'a' corresponds to the asymptote to which thresholds tend 
%   with increasing session. Parameter 'b' corresponds to the total 
%   'decay' that is, the distance between the threshold in the initial 
%   session and the asymptotic threshold (thus, a + b would correspond to 
%   the threshold at the first session). Parameter 'c' describes the decay 
%   rate or the steepness of the decline of thresholds with increasing 
%   sessions.
%
%   Thus, one might wish to define a model which constrains thresholds to 
%   follow an exponential decay function as a function of session. Rather 
%   than estimating a separate threshold for each condition, one would 
%   instead estimate the values of the three parameters of the best-fitting 
%   exponential decay function.
%
%   This requires the user to create a function which implements the
%   reparametrization. The function should accept a vector containing the 
%   parameters of the exponential decay function as its input and return
%   a vector containing the corresponding thresholds. For the situation in
%   which an observer performed in 8 consecutive sessions this function
%   could be:
%
%       function thresholds = MyThresholdReparametrization(params)
%
%       session = 1:8;
%       thresholds = params(1) + params(2)*exp(-params(3)*(session-1));
%
%   User must also create a structure (of any name) with fields '.funcA',
%   '.paramsValuesA', and '.paramsFreeA'. Field '.funcA' is assigned a 
%   handle to the reparametrization function. Field '.paramsValuesA' is 
%   assigned a vector containing initial guesses or fixed values for the 
%   parameters in the reparametrization function. Field '.paramsFreeA' is a 
%   vector of equal length to '.paramsValuesA' which contains ones for the 
%   free parameters in '.paramsFreeA' and zeros for fixed parameters. For 
%   this example, one might define:
%
%   funcParams.funcA = @MyThresholdReparametrization;
%   funcParams.paramsValuesA = [1 2 .5];
%   funcParams.paramsFreeA = [1 1 1];
%
%   Similarly, the user may define reparametrizations for the slopes, guess
%   rates and lapse rates. In order to reparametrize the slopes also, three
%   new fields need to be added to the parametrization structure: .funcB,
%   .paramsValuesB, and .paramsFreeB. For the guess rates: .funcG,
%   .paramsValuesG, and .paramsFreeG. For the lapse rates: .funcL,
%   .paramsValuesL, and .paramsFreeL. (The letters A, B, G, and L stand for
%   alpha, beta, gamma, and lambda, respectively).
%
%   To illustrate, let us also reparametrize the slopes. The following
%   reparametrization constrains the slope to be equal in all sessions and
%   also restricts this single slope estimate to be a positive number (this 
%   will avoid the fitting procedure wandering into areas of parameter 
%   space outside of the domain of the slope). We create a function that 
%   performs the reparametrization:
%
%   function slopes = MySlopeReparametrization(params)
%
%   slopes(1:8) = exp(params(1));
%
%   Palamedes will find the best-fitting value for 'params(1)' which will
%   always correspond to positive values of 'slopes'.
%
%   We also create the necessary fields in our 'funcParams' structure:
%
%   funcParams.funcB = @MySlopeReparametrization;
%   funcParams.paramsValuesB = [log(2)];
%   funcParams.paramsFreeB = [1];
%
%   Note that if we wish to custom-parametrize more than one of the PF's
%   parameters (here: thresholds and slopes) we must do so using a single 
%   reparametrization structure (here: 'funcParams').
%
%   We may now pass the structure to the above functions (instead of the
%   options 'fixed', 'constrained', 'unconstrained', or a contrast matrix).
%   For example, the call:
%
%   [paramsFitted LL exitflag output funcParamsFitted] = ...
%       PAL_PFML_FitMultiple(StimLevels,NumPos,OutOfNum, params,PF,...
%       'thresholds',funcParams,'slopes',funcParams,'guessrates',...
%       'fixed','lapserates','fixed');
%
%   constrains the thresholds and slopes according to the specifications in
%   funcParams outlined above but fixes the guess rates and lapse rates to 
%   the values specified in 'params'. Note that we may 'mix and match': 
%   for example, here we custom-parametrize thresholds and slopes but use 
%   the generic option 'fixed' for the guess rates and lapse rates. As 
%   additional examples using the 'mix and match' approach, the following 
%   calls will lead to identical fits as the above call (although these 
%   calls might lead to negative values for the slope estimate):
%   
%   [paramsFitted LL exitflag output funcParamsFitted] = ...
%       PAL_PFML_FitMultiple(StimLevels,NumPos,OutOfNum, params,PF,...
%       'thresholds',funcParams,'slopes','constrained','guessrates',...
%       'fixed','lapserates','fixed');
%
%   [paramsFitted LL exitflag output funcParamsFitted] = ...
%       PAL_PFML_FitMultiple(StimLevels,NumPos,OutOfNum, params,PF,...
%       'thresholds',funcParams,'slopes',[1 1 1 1 1 1 1 1],'guessrates',...
%       'fixed','lapserates','fixed');
%
%   The return argument 'paramsFitted' contains the threshold, slope, 
%   guess-rate, and lapse-rate estimates for each of the datasets as 
%   before. The return argument 'funcParamsFitted' is a copy of 
%   'funcParams' but with the best-fitting parameter estimates for the free 
%   parameters. The threshold estimates and slope estimates in 
%   'paramsFitted' are those that correspond to the parameter estimates in 
%   'funcParamsFitted'. For example, the slopes reported in 'paramsFitted' 
%   correspond to:
%
%   exp(funcParamsFitted.paramsValuesB)
%
%   The following example uses the above reparametrizations to fit some
%   simulated data collected using the Psi Method. This code relies on the
%   existence of the routines MyThresholdReparametrization.m and
%   MySlopeReparametrization.m as defined above.
%
%%%%%%%%%%%%%Start example code
%
%StimLevels = ...
%    [1    1.5    2    2.5    3    3.5    4;...
%     1    1.5    2    2.5    3    3.5    0;...
%     0    0.5    1    1.5    2    2.5    3;...
%     0    0.5    1    1.5    2    2.5    0;...
%     0    0.5    1    1.5    2    2.5    0;...
%     0    0.5    1    1.5    2    0      0;...
%     0    0.5    1    1.5    2    2.5    0;...
%     0.5  1      1.5  2      0    0      0];
%
%NumPos = ...
%     [0     2    33    93    13    81   102;...
%      6    52    75    14   122    53     0;...
%     25    25    67    20     7   119    59;...
%      4    46    72    16   161    21     0;...
%      9    17    96    75   122     7     0;...
%      5    80    33    80   129     0     0;...
%     13    97    23    49   126    15     0;...
%    110    31   151    24     0     0     0];
%
%OutOfNum = ...
%     [1     2    64   118    17    93   105;...
%     10    89    95    17   136    53     0;...
%     40    40    95    27     9   128    61;...
%     11    83    87    21   177    21     0;...
%     17    24   134    89   129     7     0;...
%     10   123    42    89   136     0     0;...
%     27   138    31    58   131    15     0;...
%    173    36   167    24     0     0     0];
%
%
%params(1:8,1) = 0;    %will be ignored (reparametrized by function)
%params(1:8,2) = 0;    %will be ignored (reparametrized by function)
%params(1:8,3) = .5;   %fixed value for guess rate
%params(1:8,4) = .02;  %fixed value for lapse rate
%
%funcParams.funcA = @MyThresholdReparametrization;
%funcParams.paramsValuesA = [1 2 .5];
%funcParams.paramsFreeA = [1 1 1];
%
%funcParams.funcB = @MySlopeReparametrization;
%funcParams.paramsValuesB = [log(2)];
%funcParams.paramsFreeB = [1];
%
%PF = @PAL_Logistic;
%
%[paramsFitted LL exitflag output funcParamsFitted] = ...
%  PAL_PFML_FitMultiple(StimLevels,NumPos,OutOfNum,params,PF,...
%  'thresholds',funcParams,'slopes',funcParams);
% 
% Remember that 'fixed' is the default for 'guessrates' and 'lapserates'
% and thus need not be specified explicitly.
%
%%%%%%%%%%%%%%End example code
%
%   Use of custom reparametrization in all of the functions listed at the 
%   top of this file is performed in a matter analogous to that above.
%   Those functions that perform simulations and accept the 'maxTries' and 
%   'rangeTries' arguments (i.e., PAL_PFML_BootstrapParametricMultiple,
%   PAL_PFML_BootstrapNonParametricMultiple,
%   PAL_PFML_GoodnessOfFitMultiple, and PAL_PFLR_ModelComparison) can be
%   passed a 'rangeTries' vector as before (see help comments for these
%   functions) or can be passed a 'rangeTries' structure. This structure
%   should have any or all of the fields:
%
%       .rangeTries
%       .paramsValuesA
%       .paramsValuesB
%       .paramsValuesG
%       .paramsValuesL
%
%   The format of the '.rangeTries' field should be identical to that of 
%   'rangeTries' when passed as a vector (see, for example, help comments 
%   of PAL_PFML_BootstrapParametric). The other fields should be vectors of
%   the same length as the corresponding field in the reparametrization
%   structure  (e.g., structure 'funcParams' in the 
%   example above). Each entry should specify the range of jitter to be 
%   applied to the initial guesses of parameters in the reparametrization
%   structure. Initial guesses for parameters will be selected randomly
%   from a rectangular interval centered on the values specified in the
%   reparametrization structure with range specified in the rangeTries
%   structure. For example, the example above might continue:
%
%   rangeTries.paramsValuesA = [.5 1 .2];
%   rangeTries.paramsValuesB = [1];
% 
%   [SD paramsSim LLSim converged SDfunc funcParamsSim] = ...
%      PAL_PFML_BootstrapParametricMultiple(StimLevels,OutOfNum, ...
%      paramsFitted, 400, PF,'thresholds', funcParams,'slopes',...
%       funcParams,'maxTries', 20, 'rangeTries', rangeTries);
%
%   Return arguments 'SD', 'paramsSim', 'LLSim' and 'converged' are as
%   described in help comments of PAL_PFML_BootstrapParametricMultiple and 
%   PAL_PFML_BootstrapNonParametricMultiple. 'SDfunc' is a structure 
%   containing fields '.A', '.B', '.G', and '.L' which contain the standard 
%   errors (SEs) of any parameters defined in the reparametrization 
%   structure. The field '.A'. will contain SEs for parameters defining 
%   thresholds (if any), '.B' will contains SEs for parameters defining 
%   slopes, etc. 'funcParamsSim' is a structure containing best-fitting
%   parameter estimates of parameters in the reparametrization structure
%   for each of the simulations.
%
% Introduced: Palamedes version 1.1.0 (NP)