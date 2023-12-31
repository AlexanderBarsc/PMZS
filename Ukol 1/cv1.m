clear all
close all

syms t
%% Vykresleni signalu

f = (30*pi)/(2*pi)
T = 1/f

tvek = (0:T/100:T);
w = 5*cos(30*pi*tvek+pi()/4);

figure()
plot(tvek,w);
xlabel('Cas (t)');
ylabel('w(t)');
title('Prubeh zadaneho signalu w(t)');

E = int((5*cos(30*pi*t+pi/4))^2, -inf, inf)
Wstr = (1/T)*int(5*cos(30*pi*t+pi/4), 0, T)
Pstr = (1/T)*int((5*cos(30*pi*t+pi/4))^2, 0, T)
Wef = sqrt((1/T)*int((5*cos(30*pi*t+pi/4))^2, 0, T))


%% Autokorelace 

figure()
[R,tau]= xcorr(w,'unbiased'); %unbiased pro spravnou hodnotu y
plot(T/100*tau,R) %vynasobime krokem pro spravne tau
title('Autokorelace signalu')
xlabel('\tau [s]')
ylabel('R(\tau)')

 
%% Spektralni analyza

m=-10:10;
wm=1/T*(int((5*cos(30*pi*t+pi/4)*exp(-j*m*2*pi/T*t)),t,0,T))
wwm=double(wm); %lepsi zobrazeni hodnoty
A=abs(wwm);      % 
fi=angle(wwm); 
P=A.^2;
figure()
stem(m*2*pi/T,A) %amplitudove spektrum
xlabel('\omega[rad/s]')
ylabel('|Wm|')
title('Amplitudove spektrum')
figure()
stem(m*2*pi/T,fi) %fazove spektrum
xlabel('\omega[rad/s]')
ylabel('\phi')
title('Fazove spektrum')