function G = sameconv(A, B);
%  G = sameconv(A, B);
%   
%  Causally filters A with B, giving a column vector with same height as
%  A.  (B not flipped as in standard convolution).
%
%  Convolution performed efficiently in (zero-padded) Fourier domain.


[am, an] = size(A);
[bm, bn] = size(B);
nn = am+bm-1;

G = ifft(sum(fft(A,nn).*fft(flipud(B),nn),2));
G = G(1:am,:);
