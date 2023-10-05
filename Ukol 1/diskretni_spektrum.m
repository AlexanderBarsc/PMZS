clear all;
close all;

deltat=0.001;           %vzorkovaci kmitocet
% Zadani harmonickeho signalu
T=0.5
t1=0:deltat:T;
w1=4.*cos(2*pi/T*t1-pi/8)+1;

%vykresleni harmonickeho signalu
subplot(4,3,1)
plot(t1,w1);
grid on;
xlabel('t');
ylabel('w_1(t)');
title('Harmonicky signal w_1')

% Zadani periodickeho obdelnikoveho signalu
t2a=0:deltat:0.5;
w2a(1:length(t2a))=3;
t2b=0.5+deltat:deltat:2;
w2b(1:length(t2b))=0;

t2=[t2a t2b];  %casova osa
w2=[w2a w2b];  %jedna perioda periodickeho signalu w2(t)   

%vykresleni periodickeho signalu
subplot(4,3,2)
plot(t2,w2);
grid on;
xlabel('t');
ylabel('w_2(t)');
title('Periodicky signal w_2')

% Zadani periodickeho obdelnikoveho signalu
t3a=0:deltat:1.5;
w3a(1:length(t3a))=2;
t3b=1.5+deltat:deltat:3;
w3b(1:length(t3b))=-2;

t3=[t3a t3b];  %casova osa
w3=[w3a w3b];  %jedna perioda periodickeho signalu w3(t)   

%vykresleni periodickeho signalu
subplot(4,3,3)
plot(t3,w3);
grid on;
xlabel('t');
ylabel('w_3(t)');
title('Periodicky signal w_3')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%vypocet koeficientu Fourierovy rady harmonickeho signalu 
syms k t sm %symbolicke promenne
T=0.5;      %perioda signalu

N=10;       %pocet spoctenych keficientu Fourierovy rady    
m=-N:N;   

k=2*pi/T; 	%uhlovy kmitocet

sm=1/T*(int((4.*cos(2*pi/T*t+pi/8)+1)*exp(-j*k*m*t),t,0,T)); %koeficienty Fourierovy rady
smd=double(sm);             
amp=abs(smd);           %spektrum amplitudy
faze=angle(smd);        %spektrum faze
vykon=amp.^2;           %spektrum vykonu

omega=2*pi*1/T*m;       %uhlovy kmitocet, omega=2*pi*f=2*pi/T

%vykresleni amplitudoveho, fazoveho spektra a spektra vykonu
subplot(4,3,4);
stem(omega,amp)
xlabel('\omega')
ylabel('^F^R|W_m|')
grid on;
title('Amplitudove spektrum')

subplot(4,3,5);
stem(omega,faze)
xlabel('\omega')
ylabel('\Theta_m')
grid on
title('Fazove spektrum')

subplot(4,3,6);
stem(omega,vykon)
xlabel('\omega')
ylabel('^F^R|P_m|') 
grid on
title('Spektrum vykonu')


%vypocet koeficientu Fourierovy rady periodickeho signalu 
syms k t sm %symbolicke promenne
T=2      %perioda signalu

N=10;       %pocet spoctenych keficientu Fourierovy rady    
m=-N:N;   

k=2*pi/T; 	%uhlovy kmitocet

sm=1/T*int(3*exp(-j*k*m*t),t,0,0.5); %koeficienty Fourierovy rady
smd=double(sm);             
amp=abs(smd);           %spektrum amplitudy
faze=angle(smd);        %spektrum faze
vykon=amp.^2;           %spektrum vykonu

omega=2*pi*1/T*m;       %uhlovy kmitocet, omega=2*pi*f=2*pi/T

%vykresleni amplitudoveho, fazoveho spektra a spektra vykonu
subplot(4,3,7);
stem(omega,amp)
xlabel('\omega')
ylabel('^F^R|W_m|')
grid on;
title('Amplitudove spektrum')

subplot(4,3,8);
stem(omega,faze)
xlabel('\omega')
ylabel('\Theta_m')
grid on
title('Fazove spektrum')

subplot(4,3,9);
stem(omega,vykon)
xlabel('\omega')
ylabel('^F^R|P_m|') 
grid on
title('Spektrum vykonu')

%vypocet koeficientu Fourierovy rady periodickeho signalu 
syms k t sm %symbolicke promenne
T=3      %perioda signalu

N=10;       %pocet spoctenych keficientu Fourierovy rady    
m=-N:N;   

k=2*pi/T; 	%uhlovy kmitocet

sm=1/T*(int(2*exp(-j*k*m*t),t,0,1.5)+int(-2*exp(-j*k*m*t),t,1.5,3)); %koeficienty Fourierovy rady
smd=double(sm);             
amp=abs(smd);           %spektrum amplitudy
faze=angle(smd);        %spektrum faze
vykon=amp.^2;           %spektrum vykonu

omega=2*pi*1/T*m;       %uhlovy kmitocet, omega=2*pi*f=2*pi/T

%vykresleni amplitudoveho, fazoveho spektra a spektra vykonu
subplot(4,3,10);
stem(omega,amp)
xlabel('\omega')
ylabel('^F^R|W_m|')
grid on;
title('Amplitudove spektrum')

subplot(4,3,11);
stem(omega,faze)
xlabel('\omega')
ylabel('\Theta_m')
grid on
title('Fazove spektrum')

subplot(4,3,12);
stem(omega,vykon)
xlabel('\omega')
ylabel('^F^R|P_m|') 
grid on
title('Spektrum vykonu')
