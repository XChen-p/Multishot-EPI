function M = fftdim(M,dim)

%
% [m] = fftdim(M,dim)
%
% performs centric Fourier Transform along dimension dim.
%

for i = dim
    M   =   fftshift(fft(ifftshift(M, i), [], i), i)/sqrt(size(M,i));
end