clear all
close all

%% Cviceni 3 vypracovani

% Vykreslene signaly m1 a m2

tsampling = 0.001

t = 0:tsampling:1
f1 = 5;
omega1 = 2*pi*f1;
U1 = 2;

f2 = 8;
omega2 = 2*pi*f2;
U2 = 3;

m1 = U1*cos(omega1*t);

m2 = U2*cos(omega2*t);

figure();
plot(t,m1);
hold on
plot(t,m2);
ylabel('Amplituda');
xlabel('Cas (t)');
title('Vykreslene signaly m1 a m2');

%% Soucet signalu m1 a m2 -> m3

m3 = m1 + m2;

tmaxharmonic = 1/f2;
tsamplingNew = tmaxharmonic/5;

tfit = tsamplingNew / tsampling

t2periods = (2*tmaxharmonic)/tsampling


newValues = m3(1:tfit:t2periods); % Vytahneme kazdy petinovy vzorek z periody