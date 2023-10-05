%Machacek

%Korelacni funkce

clear all;
close all;

%definice sinusoveho signalu w1(t)

deltat=0.001; %krok

t=0:deltat:4-deltat;            %cas
w1=5.*sin(4*pi*t+pi/8)+2;     %funkce y(t)

%vykresleni funkce w1(t) do grafickeho okna c.1
figure(1),
subplot(351);
plot(t,w1)     
xlabel('t [s]')         %popis x-ove osy
ylabel('w_1(t)')          %popis y-ove osy
title('Sinusovy signal')

%Definice obdelníkoveho signalu w2(t)
t2a=0:deltat:2-deltat;
w2a(1:length(t2a))=0; %signal w1(t) na intervalu <0,1>

t2b=2:deltat:3;
w2b(1:length(t2b))=2; %signal w1(t) na intervalu <1,3>

t2c=3+deltat:deltat:4;
w2c(1:length(t2c))=0; %signal w1(t) na intervalu <3,4>

t2=[t2a t2b t2c]; % casova osa
w2=[w2a w2b w2c]; % signal w1(t)

%vykresleni funkce w2(t)
subplot(352);
plot(t2,w2);

xlabel('t[s]')
ylabel('w_2(t)')
title('Obdelnikovy signal')

%definice cosinusoveho signalu w3(t)
w3=2.*cos(pi*t)-1;     %funkce y(t)

%vykresleni funkce w3(t) do grafickeho okna c.3
subplot(353);
plot(t,w3)     
xlabel('t [s]')         %popis x-ove osy
ylabel('w_3(t)')          %popis y-ove osy
title('Kosinusovy signal')

%definice nahodneho signalu w4(t)
w4=10*randn(1,length(t));     %funkce y(t)

%vykresleni funkce w3(t) do grafickeho okna c.3
subplot(354);
plot(t,w4)     
xlabel('t [s]')         %popis x-ove osy
ylabel('w_4(t)')          %popis y-ove osy
title('Nahodny signal')


%Definice periodickeho obdelníkoveho signalu w41(t)
sirka=0.25;                                      %šíøka pulsu
t41=0;                                           %????
w41=0;                                           %????
perioda =0.5;                                     %perioda signálu    
delka=4;   
for n=0:perioda:delka-perioda+2*deltat;
    t41a=n:deltat:n+sirka;
    w41a(1:length(t41a))=1;
    t41b=n+sirka+deltat:deltat:n+perioda-deltat;
    w41b(1:length(t41b))=0;

    t41=[t41 t41a t41b]; % casova osa 
    w41=[w41 w41a w41b]; % signal p(t)

end;

%vykresleni funkce w41(t)
subplot(355);
plot(t41,w41);

xlabel('t[s]')
ylabel('w_4_1(t)')
title('Periodicky Obdelnikovy signal')

%energeticky sinusovy signal w42(t)
t42a=0:deltat:1-deltat;
w42a(1:length(t42a))=0; %signal w6(t) na intervalu <0,1>

t42b=1:deltat:3;
w42b(1:length(t42b))=10.*sin(2*pi*t42b+pi/8)+3; %signal w42(t) na intervalu <1,3>

t42c=3+deltat:deltat:4;
w42c(1:length(t42c))=0; %signal w42(t) na intervalu <3,4>

t42=[t42a t42b t42c]; % casova osa
w42=[w42a w42b w42c]; % signal w42(t)

%vykresleni funkce w6(t)
subplot(356);
plot(t42,w42);

xlabel('t[s]')
ylabel('w_4_2(t)')
title('energeticky sinusovy signal')

%Soucet harmonickych signalu
w5=w1+w3;
subplot(357);
plot(t,w5);

xlabel('t[s]')
ylabel('w5(t)')
title('Soucet harmon. signalu')

%Soucet harmonickeho signalu a sumu
w6=w1+w4;
subplot(358);
plot(t,w6);

xlabel('t[s]')
ylabel('w6(t)')
title('Soucet harmon.+sum signalu')


%Vypocet autokorelacni funkce energetickeho signalu
Ra1=deltat*xcorr(w2,w2); 
subplot(359);
plot(-deltat*(length(Ra1)-1)/2:deltat:deltat*(length(Ra1)-1)/2,Ra1);

xlabel('tau[s]')
ylabel('Ra(tau)')
title('Autok. energ. obd. signalu')

Ra2=deltat*xcorr(w42,w42); 
subplot(3,5,10);
plot(-deltat*(length(Ra2)-1)/2:deltat:deltat*(length(Ra2)-1)/2,Ra2);

xlabel('tau[s]')
ylabel('Ra(tau)')
title('Autok. energ. sin. signalu')

%Vypocet autokorelacni funkce vykonoveho signalu
[Ra3,tau1]=xcorr(w1,length(w1),'unbiased'); 
tau=tau1*deltat;
subplot(3,5,11);
plot(tau,Ra3);

xlabel('tau[s]')
ylabel('Ra(tau)')
title('Autok. vykon. harm. signalu')

%Vypocet autokorelacni funkce vykonoveho signalu
[Ra5,tau2]=xcorr(w5,length(w5),'unbiased'); 
tau3=tau2*deltat;
subplot(3,5,12);
plot(tau3,Ra5);

xlabel('tau[s]')
ylabel('Ra(tau)')
title('Autok. souctu vykon. harm. signalu')

%Vypocet autokorelacni funkce vykonoveho signalu
[Ra6,tau4]=xcorr(w6,length(w6),'unbiased'); 
tau5=tau4*deltat;
subplot(3,5,13);
plot(tau5,Ra6);

xlabel('tau[s]')
ylabel('Ra(tau)')
title('Autok. souctu vykon. harm.+sum')

%Vypocet vzajemna korelacni funkce vykonovych signalu w1 a w3
[Rab7,tau1]=xcorr(w1,w3,length(w1),'unbiased'); 
tau=tau1*deltat;
subplot(3,5,14);
plot(tau,Rab7);

xlabel('tau[s]')
ylabel('Rab(tau)')
title('Vzajemna kor. fce vykon. signalu')

%Vypocet vzajemna korelacni funkce energetickych signalu w2 a w42
[Rab8,tau1]=xcorr(w2,w42,length(w1),'biased'); 
tau=tau1*deltat;
subplot(3,5,15);
plot(tau,Rab8);

xlabel('tau[s]')
ylabel('Rab(tau)')
title('Vzajemna kor. fce energ. signalu')