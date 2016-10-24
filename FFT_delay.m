function [ sd ] = FFT_delay(s,d)
% this funtion accepts 2 arguments:
% 1. input signal s
% 2. time delay in seconds

% return the delay version of the signal sd(column vector)

S=fftshift(fft(s)); % get the frequency content of the signal
N=length(s);
%omega=(-pi+2*pi/N):(2*pi/N):pi;
omega = -pi:(2*pi/N):(pi-2*pi/N); % define the frequency bins
for p=1:N
    S(p)=S(p)*exp(-1i*omega(p)*d); % multiply each frequency bin with appropriate weights
end

sd=real(ifft(ifftshift(S)));
sd = sd';
end






