%Priklad 1
clear all;
close all;

%Vykresleni funkce y(t)
t=-5:0.01:5;            %cas
y=5.*cos(2*t+pi/4);     %funkce y(t)
figure(1),plot(t,y)     %vykresleni funkce y(t) do grafickeho okna c.1
xlabel('t [s]')         %popis x-ove osy
ylabel('y(t)')          %popis y-ove osy

%Zmena rozsahu y-ove osy
%figure(2),axis([0 pi -10 10])   %zmena rozsahu y-ove souradnice

%Vypocet funkce r(t)=y(t)*u(t)
u=sin(4*t)+8;
r=y.*u;                 %nasobeni dvou vektoru, prvek po prvku
figure(3),subplot(3,1,1),plot(t,y,'b')  %v grafickem okne c.3 jsou zobrazeny 
subplot(3,1,2),plot(t,u,'g')            %3 ruzne grafy
subplot(3,1,3),plot(t,r,'r')

%Vypocet komplexni funkce
C=sqrt(2)*(1+j);
a=1+2*pi*j;
y=C.*exp(a*t);
re_y=real(y);           %realna cast komplexniho cisla
im_y=imag(y);           %imaginarni cast komplexniho cisla
figure(4),plot(re_y,im_y)  

%Ulozeni promennych do souboru cviceni_2.mat
save cviceni_2 re_y im_y    %ulozeni promennych re_y, im_y 
                            



%Priklad 2

t=0:0.01:8;                 		%casova promenna
y1=2.*sin(0.5*pi*t+pi/2);   		%zadany signal
figure(5), plot(t,y1)                  		%vykresleni signalu
xlabel('t [s]')             		%popis x-ove souradnice
ylabel('y1(t)')             		%popis y-ove souradnice    
title('y1(t)=2*sin(0.5*pi*t+pi/2)')	%titulek grafu




%Priklad 3

T=0.25;		% perioda signalu ya(t) je T=0.25

syms t		% zavedeni symbolicke promenne t
str=1/T*int(3*cos(8*pi*t+pi/4),0,T) 	%stredni hodnota

