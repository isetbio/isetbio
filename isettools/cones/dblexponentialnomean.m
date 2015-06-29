function [fit] = dblexponentialnomean(coef, x)

fit = coef(1) .* exp(-x * abs(coef(2))) + coef(3) .* exp(-x * abs(coef(4)));
