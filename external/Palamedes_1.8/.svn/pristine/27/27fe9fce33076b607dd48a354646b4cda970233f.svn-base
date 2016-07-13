% ParametrizeThresholds Reparametrizes thresholds of PFs for
%   PAL_PFLR_CustomDefine_Demo.
%
%NP (November 2009)

function alphas = ParametrizeThresholds(params)

alphas(1:8) = (params(1)-params(4)/2) + ...
  (params(2)-params(5)/2)*exp(-(params(3)-params(6)/2)*[0:7]);
alphas(9:16) = (params(1)+params(4)/2) + ...
  (params(2)+params(5)/2)*exp(-(params(3)+params(6)/2)*[0:7]);