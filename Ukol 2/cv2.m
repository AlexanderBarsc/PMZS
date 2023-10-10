clear all
close all

%% Vykresleni diskretniho signalu

N = 10; % Pocet vzorku
Fs = 25; % Vzorkovaci frekvence
Ts = 1/Fs; % Vzorkovaci perioda

Wmax = 0.2; % Amplituda
omega = 10*pi; % Uhlovy kmitocet
phase = -pi; % Fazovy posun

n = 0:1:(N - 1); % N potrebne pro jednu periodu

w = Wmax*cos(omega*n*Ts+pi) + 0.5
stem(n,w);
xlabel('n');
ylabel('w[n]');
title('Harmonický signál w[n]');


W_k1=fft(w);
W_k=fftshift(W_k1);

k = -(N/2):1:N/2-1;

figure;
subplot(2,1,1)
stem(k,(abs(W_k)/N))
xlabel('k[-]')
ylabel('|W[k]|')
title(' Amplitude frequency spectra |W[k]|')

subplot(2,1,2)
stem(k,(atan2(imag(W_k),real(W_k))))
xlabel('k[-]')
ylabel('phase W[k]')
title(' Phase frequency spectra W[k]')

%% Obdelnikova okenni funkce
