function ts_shuffled = phase_shuffling(ts)
% Add a random phase to single time serie
%
% Inputs
%    ts - time serie (1D)
%
% Outputs
%    ts_shuffled - surrogate time serie
%

if length(ts) ~= numel(ts); error('Input must to be single time serie'); end
dim = size(ts);
ts = ts(:);
time = length(ts);
if mod(time, 2) == 0; time = time - 1; ts = ts(1:time); end

% Fourier transform for time serie
fft_ts = fft(ts);

% Random phase
random_phase = rand(floor(length(fft_ts)/2), 1);

% Fourier anti-transform for real time serie
random_phase_function = [exp(2*pi*1i*random_phase); conj(flipud(exp(2*pi*1i*random_phase)))];
fft_ts_shuffled = fft_ts;
fft_ts_shuffled(2:end) = fft_ts(2:end).*random_phase_function;
ts_shuffled = real(ifft(fft_ts_shuffled));
if dim(1) < dim(2); ts_shuffled = ts_shuffled'; end

end