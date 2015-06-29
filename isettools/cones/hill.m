function [fit] = hill(beta, x)

%fit = log10(1.0 ./ (1 + (beta(1) ./ x).^beta(2)));
fit = 1.0 ./ (1 + (beta(1) ./ x).^beta(2));
