function [fit] = weber_fechner(wf, x)

fit = log10(1 ./ (1 + x ./ wf(1)));
