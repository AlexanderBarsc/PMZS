
clear all;
close all;
clc;

%--------------------------------------------------------------------------
%Example 3_1
%--------------------------------------------------------------------------
clear all;
close all;

%--------------------------------------------------------------------------
%Signal definition w[n]
N=10;            %signal period 
POINT_NUMBER = N; %period numbers
Ts=1/25; %sampling frequency
% Harmonic signal defintion
Wmax = 2;
OMEGA = 10*pi;
phase = -pi;

n=0:1:POINT_NUMBER-1;				%time axe
w=Wmax*cos(OMEGA*n*Ts+phase)+1;		% signal w[n]

%% signal

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




%% fft



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



%% window 1

wo = ones(1,POINT_NUMBER);
wa = w.*wo;
figure;
%Signal graphs w[n]
subplot(2,1,1)
stem(n,wa); 
xlabel('n')
ylabel('wa[n]')
title('Signal wa[n]')

%Signal graphs wo[n]
t1 = Ts.*n;
subplot(2,1,2)
stem(t1,wa); 
xlabel('t[s]')
ylabel('wa[t]')
title(['Signal wa[t]'])

%% fft with window 1


wa = w.*wo;


Wa_k1=fft(wa);
Wa_k=fftshift(Wa_k1);

k = -(POINT_NUMBER/2):1:POINT_NUMBER/2-1;

figure;
subplot(2,1,1)
stem(k,(abs(Wa_k)/POINT_NUMBER))
xlabel('k[-]')
ylabel('|Wa[k]|')
title(' Amplitude frequency spectra |Wa[k]|')

subplot(2,1,2)
stem(k,(atan2(imag(Wa_k),real(Wa_k))))
xlabel('k[-]')
ylabel('phase Wa[k]')
title(' Phase frequency spectra Wa[k]')

%% window 2




wo2=1-cos(2*pi*n/(N-1));	


wa2 = w.*wo2;
figure;
%Signal graphs w[n]
subplot(2,1,1)
stem(n,wa2); 
xlabel('n')
ylabel('wa2[n]')
title('Signal wa2[n]')

%Signal graphs wo2[n]
t1 = Ts.*n;
subplot(2,1,2)
stem(t1,wa2); 
xlabel('t[s]')
ylabel('wa2[t]')
title(['Signal wa2[t]'])


%% fft window 2



Wa2_k1=fft(wa2);
Wa2_k=fftshift(Wa2_k1);

k = -(POINT_NUMBER/2):1:POINT_NUMBER/2-1;

figure;
subplot(2,1,1)
stem(k,(abs(Wa2_k)/POINT_NUMBER))
xlabel('k[-]')
ylabel('|Wa2[k]|')
title(' Amplitude frequency spectra |Wa2[k]|')

subplot(2,1,2)
stem(k,(atan2(imag(Wa2_k),real(Wa2_k))))
xlabel('k[-]')
ylabel('phase Wa2[k]')
title(' Phase frequency spectra Wa2[k]')


%% window 3

wo3=1-1.98*cos(2*pi*n/(N-1))+1.29*cos(4*pi*n/(N-1))-0.388*cos(6*pi*n/(N-1))+0.0322*cos(8*pi*n/(N-1));	



wa3 = w.*wo3;
figure;
%Signal graphs w[n]
subplot(2,1,1)
stem(n,wa3); 
xlabel('n')
ylabel('wa3[n]')
title('Signal wa3[n]')

%Signal graphs wo3[n]
t1 = Ts.*n;
subplot(2,1,2)
stem(t1,wa3); 
xlabel('t[s]')
ylabel('wa3[t]')
title(['Signal wa3[t]'])

%% fft 3




Wa3_k1=fft(wa3);
Wa3_k=fftshift(Wa3_k1);

k = -(POINT_NUMBER/2):1:POINT_NUMBER/2-1;

figure;
subplot(2,1,1)
stem(k,(abs(Wa3_k)/POINT_NUMBER))
xlabel('k[-]')
ylabel('|Wa3[k]|')
title(' Amplitude frequency spectra |Wa3[k]|')

subplot(2,1,2)
stem(k,(atan2(imag(Wa3_k),real(Wa3_k))))
xlabel('k[-]')
ylabel('phase Wa3[k]')
title(' Phase frequency spectra Wa3[k]')


%% window 4



n1 = 0:1:(POINT_NUMBER-1)/2;
n2 = (POINT_NUMBER-1)/2 +1/2: 1 : POINT_NUMBER-1;



wo4_1 = 2*n1/(POINT_NUMBER-1);
wo4_2 = 2-2*n2/(POINT_NUMBER-1);
wo4 = [wo4_1, wo4_2];


wa4 = w.*wo4;
figure;
%Signal graphs w[n]
subplot(2,1,1)
stem(n,wa4); 
xlabel('n')
ylabel('wa4[n]')
title('Signal wa4[n]')

%Signal graphs wo4[n]
t1 = Ts.*n;
subplot(2,1,2)
stem(t1,wa4); 
xlabel('t[s]')
ylabel('wa4[t]')
title(['Signal wa4[t]'])


%% fft 4




Wa4_k1=fft(wa4);
Wa4_k=fftshift(Wa4_k1);

k = -(POINT_NUMBER/2):1:POINT_NUMBER/2-1;

figure;
subplot(2,1,1)
stem(k,(abs(Wa4_k)/POINT_NUMBER))
xlabel('k[-]')
ylabel('|Wa4[k]|')
title(' Amplitude frequency spectra |Wa4[k]|')

subplot(2,1,2)
stem(k,(atan2(imag(Wa4_k),real(Wa4_k))))
xlabel('k[-]')
ylabel('phase Wa4[k]')
title(' Phase frequency spectra Wa4[k]')


%% window 5

wo5 = 0.5 - 0.5*cos(2*pi*n/(N-1));

wa5 = w.*wo5;
figure;
%Signal graphs w[n]
subplot(2,1,1)
stem(n,wa5); 
xlabel('n')
ylabel('wa5[n]')
title('Signal wa5[n]')

%Signal graphs wo5[n]
t1 = Ts.*n;
subplot(2,1,2)
stem(t1,wa5); 
xlabel('t[s]')
ylabel('wa5[t]')
title(['Signal wa5[t]'])


%% fft 5




Wa5_k1=fft(wa5);
Wa5_k=fftshift(Wa5_k1);

k = -(POINT_NUMBER/2):1:POINT_NUMBER/2-1;

figure;
subplot(2,1,1)
stem(k,(abs(Wa5_k)/POINT_NUMBER))
xlabel('k[-]')
ylabel('|Wa5[k]|')
title(' Amplitude frequency spectra |Wa5[k]|')

subplot(2,1,2)
stem(k,(atan2(imag(Wa5_k),real(Wa5_k))))
xlabel('k[-]')
ylabel('phase Wa5[k]')
title(' Phase frequency spectra Wa5[k]')










