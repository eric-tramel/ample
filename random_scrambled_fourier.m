function A = random_scrambled_fourier(N,M,varargin)
% RANDOM_SCRAMBLED_FOURIER Generate a random projection operator by
% scrambling and sub-sampling of the fourier space.
scram = randperm(N);
[~,inv_scram] = sort(scram);
I = randperm(N);
I = I(1:M);

% Handle user inputs
while ~isempty(varargin)
    switch lower(varargin{1})
        case 'scram_pattern'
            scram = varargin{2};
            [~,inv_scram] = sort(scram);
        case 'row_select'
            I = varargin{2};
    end    
    varargin(1:2) = [];
end

% Set the operator values
A.N = N;
A.M = M;
A.forward = @(x_) rsfft_forward(x_,I,scram);
A.adjoint = @(x_) rsfft_adjoint(x_,N,I,inv_scram);
A.squared_forward = @(x_) rsfft_squared_forward(x_,M);
A.squared_adjoint = @(x_) rsfft_squared_adjoint(x_,N);

function y = rsfft_forward(x,I,scram)
% RANDOM_SCRAMBED_FOURIER::RSFFT_FORWARD Perform the forward projection
% operation.
X = fft(x(scram));              % Get the fourier coefficients of the scrambled signal
y = X(I);                       % Sub-sample the spectra (i.e. M random rows of FFT).


function x = rsfft_adjoint(y,N,I,inv_scram)
% RANDOM_SCRAMBED_FOURIER::RSFFT_ADJOINT Perform the adjoint operation via
% the IFFT.
Y    = zeros(N,1);
Y(I) = y;                     % Inserting the Y's into place  
x = ifft(Y);                  % Perform the IFFT
x = x(inv_scram);             % Reverse the scrambling
x = real(x);                  % Should we make this real? 

function y = rsfft_squared_forward(x,M)
% RANDOM_SCRAMBLED_FOURIER::RSFFT_SQUARED_FORWARD Perform the
% squared-forwad operator. Since this is a Fourier operation, the squaring
% of all the terms (which includes an absolute-value computation) forces
% all entries in the 'explicit' Fourier transform matrix to go to 1. This
% means that all of the M entries are now the sum of x.
y = mean(x).*ones(M,1);     % Actually, these should be means


function x = rsfft_squared_adjoint(y,N)
% RANDOM_SCRAMBLED_FOURIER::RSFFT_SQUARED_FORWARD Perform the
% squared-forwad operator. 
x = mean(y).*ones(N,1);    % Actually, these should be means
