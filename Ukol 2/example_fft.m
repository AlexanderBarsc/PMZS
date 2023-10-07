%--------------------------------------------------------------------------
%Example Cv.2 - fft
%--------------------------------------------------------------------------
clear all;
close all;

%--------------------------------------------------------------------------
%Signal definition w[n]
N=10;            %signal period 
POINT_NUMBER = N; %period numbers
Ts=1/50; %sampling frequency
% Harmonic signal defintion
Wmax = 4;
OMEGA = 10*pi;
phase = pi;

n=0:1:POINT_NUMBER-1;				%time axe
w=Wmax*cos(OMEGA*n*Ts+phase);		% signal w[n]

%Signal graphs w[n]
subplot(2,1,1)
stem(n,w); 
xlabel('n')
ylabel('w[n]')
title('Harmonic signal w[n]')

%Signal graphs w[n]
t1 = Ts.*n;
subplot(2,1,2)
stem(t1,w); 
xlabel('t[s]')
ylabel('w[t]')
title(' Harmonic signal w[t]')

W_k1=fft(w);
W_k=fftshift(W_k1);

k = -(POINT_NUMBER/2):1:POINT_NUMBER/2-1;

figure;
subplot(2,1,1)
stem(k,(abs(W_k)/POINT_NUMBER))
xlabel('k[-]')
ylabel('|W[k]|')
title(' Amplitude frequency spectra |W[k]|')

subplot(2,1,2)
stem(k,(atan2(imag(W_k),real(W_k))))
xlabel('k[-]')
ylabel('phase W[k]')
title(' Phase frequency spectra W[k]')

