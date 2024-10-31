% understand_fft
clear all;
clf('reset');

% option A to set time and frequency bins: nice round numbers
dt = 10^(-10); % 0.1 ns, time between samples
fs = 1/dt; % sampling freq, 10 GHz
N = 1000; % number of time samples
 
% option B: match the VNA
fdiff = (10^7)/3; % svna from 700 MHz - 3 GHz w/ 750 data points
N = 2000; % N > 2*750
fs = N*fdiff;
dt = 1/fs;
t = 0:dt:(N-1)*dt;

% option 1: cos, sin, and noise to play with
freq = 1*10^9;
t = 0:dt:(N-1)*dt;
c = cos(2*pi*freq*t);
s = 0.5 * sin(2*pi*(2*freq)*(t-33*dt));
n = randn(1,N) * 100;
a = c + s + n;

% option 2: simulated radar pulse returns
a = zeros(1,N);
a(1) = 1;
%a = a + randn(1,N) * .001;

%plot(t,a);
%title('TD a');

A = fft(a); % complex amplitudes of frequency spectrum of A, native ordering
Ashift = fftshift(A); % orders the frequencies to match f
fshift = -fs/2:fs/N:fs/2;
fshift = fshift(1:N);
f = fftshift(fshift); % this is ordered in the "native" order

%figure
%plot(fshift,abs(Ashift));
%hold on;
%plot(fshift,angle(Ashift));
%title('A shift amplitude and phase')
%figure 
%plot(fshift,real(Ashift));
%hold on
%plot(fshift,imag(Ashift));
%title('A shift real and imaginary');
%hold off;

Bshift = Ashift;
% option I: do nothing
% option II: modify B by masking out certain frequencies
%mask = abs(fshift) < 1.5*10^9;
%mask = abs(fshift) < 0.1*10^9;
%mask = (abs(fshift) > 700*10^6) .* (abs(fshift) < 3*10^9);
mask = (abs(fshift) > 700*10^6) .* (abs(fshift) < 3*10^9);
Bshift = mask.*Bshift;

% option III: overwrite Bshift with VNA data in the correct slots
% Bshift = 0*Bshift;
% data = svna_data_analysis(15);
% vna_start = 1152; % index of 503.33 MHz
% vna_end = 1901; % index of 3 GHz
% vna_negative_start = 101; % -3
% vna_negative_end = 850; % -503.3
% vna_comp = data(4,:) + j*data(5,:);
% Bshift(vna_start:vna_end) = vna_comp;  % insert the 750 real and imag samples into the positive frequencies
% Bshift(vna_negative_start:vna_negative_end) = flip(conj(vna_comp));  % do the same to negative freqs.

%figure
plot(fshift,real(Bshift));
hold on
plot(fshift,imag(Bshift));
title('B shift real and imaginary');

B = fftshift(Bshift);
b = ifft(B);
figure 
plot(t,real(b));
hold on
plot(t,imag(b));
title('TD b real and imaginary');

figure
plot(fshift,B);
hold on
plot(fshift,Bshift);
