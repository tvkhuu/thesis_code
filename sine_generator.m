function [x] = sine_generator(Amplitude,num_of_sam,norm_freq)

% this function generates a sine wave with the amplitude,number of samples,

% frequency and sampling frequency.

% assume no phase angle

n = 0:num_of_sam-1;

freq = 500; %Hz

upsample_factor = 1/norm_freq;

fs = (2*freq)*upsample_factor;

x=Amplitude*sin(2*pi*freq*n/fs);

%{

% for plotting

tScaleinSecond = n/fs;   

plot(tScaleinSecond, x)

xlabel('time in second') 

ylabel('magnitude')

title('sine wave')

%}

end