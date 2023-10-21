%--------------------------------------------------------------------------
% Vy�i�t�n� v�ech prom�nn�ch pam�ti a uzav�en� v�ech grafick�ch oken
clear all;
close all;

%--------------------------------------------------------------------------
% Specifikace �asov�ho pr�b�hu-�as t[s] sign�lu
cas = 10; % �as pr�b�hu simulace 
pocet_bodu = 100; %po�et bod� v �ase b�hem jedn� vzorkovac� periody fs
T_s = 0.5; % Perioda vzorkov�n�
deltat = T_s/pocet_bodu; % zanedbateln� �asov� okam�ik
t = 0:deltat:cas; % �asov� pr�b�h
%--------------------------------------------------------------------------
% 1. zad�n� harmonick�ho signalu w1(t)
W_max = 10; % Amplituda sign�lu
f = 1/4; % frekvence sign�lu
omega = 2*pi*f; % �hlov� kmito�et 
fi = pi; % f�zov� posun
C = 1; % konstantn� slo�ka
w1 = W_max*cos(omega*t+fi) + C; % harmonick� signal w1(t) 
%--------------------------------------------------------------------------
% 2. zad�n� obdeln�kov�ho signalu w2(t)
t2a=0:deltat:2-deltat;
w2a(1:length(t2a))=0; %1.��st signalu w2(t)

t2b=2:deltat:4;
w2b(1:length(t2b))=3; %2.��st signalu w2(t)

t2c=4+deltat:deltat:cas;
w2c(1:length(t2c))=0; %3.��st signalu w2(t)

w2=[w2a w2b w2c]; % slo�en� celkov� signal w2(t)
%--------------------------------------------------------------------------
% 3. zad�n� konstantn�ho sign�lu w3(t)
w3(1:length(t))=5; %konstantn� signal w3(t) 

%--------------------------------------------------------------------------
% V�PO�ET IMPULSN�HO SIGN�LU - IDE�LN� VZORKOV�N�
%--------------------------------------------------------------------------
% V�po�et posloupnosti Diracov�ch impulsu
t4=0;  
w4=0;

for n=0:T_s:cas-T_s+2*deltat;
t4a=n:deltat:n;
w4a(1:length(t4a))=1;

t4b=n+deltat:deltat:n+T_s-deltat;
w4b(1:length(t4b))=0;

w4=[w4 w4a w4b]; % signal delta(t)
end;

w1_i=w1.*w4; %PAM - ide�ln� vzorkov�n� harmonick�ho sign�lu
w2_i=w2.*w4; %PAM - ide�ln� vzorkov�n� obdeln�kov�ho sign�lu
w3_i=w3.*w4; %PAM - ide�ln� vzorkov�n� konstantn�ho sign�lu

%--------------------------------------------------------------------------
% V�PO�ET IMPULSN�HO SIGN�LU - P�IROZEN� VZORKOV�N�
%--------------------------------------------------------------------------
%Definov�n� umpuls� dan� ���ky
sirka = 0.2;  % ���ka impulsu
t5=0;
w5=0;

for n=0:T_s:cas-T_s+2*deltat;
t5a=n:deltat:n+sirka;
w5a(1:length(t5a))=1;

t5b=n+sirka+deltat:deltat:n+T_s-deltat;
w5b(1:length(t5b))=0;

w5=[w5 w5a w5b]; % posloupnost impulsn�ch obdeln�kov�ch signalu p(t)
end;

w1_p=w1.*w5; %PAM - p�irozen� vzorkov�n� harmonick�ho sign�lu
w2_p=w2.*w5; %PAM - p�irozen� vzorkov�n� obdeln�kov�ho sign�lu
w3_p=w3.*w5; %PAM - p�irozen� vzorkov�n� konstantn�ho sign�lu

%--------------------------------------------------------------------------
% V�PO�ET IMPULSN�HO SIGN�LU - OKAM�IT� VZORKOV�N�
%--------------------------------------------------------------------------
%Definov�n� jednoho impulsu dan� ���ky
t6a=0:deltat:sirka;
w6a(1:length(t6a))=1;

t6b=sirka+deltat:deltat:cas;
w6b(1:length(t6b))=0;

w6=[w6a w6b]; % impulsn� obdeln�kov� signal p1(t)

w1_o1=conv(w1_i,w6); %PAM - okam�it� vzorkov�n� harmonick�ho sign�lu
w1_o = w1_o1(1:length(t));
w2_o1=conv(w2_i,w6); %PAM - okam�it� vzorkov�n� obdeln�kov�ho sign�lu
w2_o = w2_o1(1:length(t));
w3_o1=conv(w3_i,w6); %PAM - okam�it� vzorkov�n� konstantn�ho sign�lu
w3_o = w3_o1(1:length(t));

%--------------------------------------------------------------------------
% GRAFY - ZOBRAZEN� PAM MODULACE HARMONICK�HO SIGN�LU V �ASOV� OBLASTI
%--------------------------------------------------------------------------
% Vykreslen� grafick�ho �asov�ho pr�b�hu harmonick�ho sign�lu
figure;
subplot(4,1,1)
plot(t,w1); 
title('Harmonick� sign�l w_1(t)');
ylabel('w_1(t)')
xlabel('t');
grid on;

%Vykreslen� grafick�ho �asov�ho pr�b�hu posloupnosti Diracov�ch impulsu
subplot(4,2,3);
plot(t,w4);
title('Posloupnost Diracov�ch impuls�');
ylabel('\delta(t)');
xlabel('t');
grid on;

%Vykreslen� grafick�ho �asov�ho pr�b�hu impulsn�ho sign�lu vytvo�en�ho
%ide�ln�m vzorkov�n�m
subplot(4,2,4);
plot(t,w1_i);
title('PAM modulace - ide�ln� vzorkov�n�');
ylabel('w_1_i(t)')
xlabel('t');
grid on;

%Vykreslen� grafick�ho �asov�ho pr�b�hu posloupnosti obdeln�kov�ch impuls�
subplot(4,2,5);
plot(t,w5);
title('Posloupnost posloupnosti obdeln�kov�ch impulsu');
ylabel('p(t)');
xlabel('t');
grid on;

%Vykreslen� grafick�ho �asov�ho pr�b�hu impulsn�ho sign�lu vytvo�en�ho
%p�irozen�m vzorkov�n�m
subplot(4,2,6);
plot(t,w1_p);
title('PAM modulace - p�irozen� vzorkov�n�');
ylabel('w_1_i_p(t)')
xlabel('t');
grid on;

%Vykreslen� grafick�ho �asov�ho pr�b�hu obdeln�kov�ho impulsu
subplot(4,2,7);
plot(t,w6);
title('Obdeln�kov� impuls');
ylabel('p_1(t)');
xlabel('t');
grid on;

%Vykreslen� grafick�ho �asov�ho pr�b�hu impulsn�ho sign�lu vytvo�en�ho
%okam�it�m vzorkov�n�m
subplot(4,2,8);
plot(t,w1_o);
title('PAM modulace - okam�it� vzorkov�n�');
ylabel('w_1_i_o(t)')
xlabel('t');
grid on;

%--------------------------------------------------------------------------
% GRAFY - ZOBRAZEN� PAM MODULACE PERIODICK�HO OBDELN�KOV�HO SIGN�LU V �ASOV� OBLASTI
%--------------------------------------------------------------------------
% Vykreslen� grafick�ho �asov�ho pr�b�hu periodick�ho obdeln�kov�ho sign�lu
figure;
subplot(4,1,1)
plot(t,w2); 
title('Obdeln�kov� periodick� sign�l w_2(t)');
ylabel('w_2(t)')
xlabel('t');
grid on;

%Vykreslen� grafick�ho �asov�ho pr�b�hu posloupnosti Diracov�ch impulsu
subplot(4,2,3);
plot(t,w4);
title('Posloupnost Diracov�ch impuls�');
ylabel('\delta(t)');
xlabel('t');
grid on;

%Vykreslen� grafick�ho �asov�ho pr�b�hu impulsn�ho sign�lu vytvo�en�ho
%ide�ln�m vzorkov�n�m
subplot(4,2,4);
plot(t,w2_i);
title('PAM modulace - ide�ln� vzorkov�n�');
ylabel('w_2_i(t)')
xlabel('t');
grid on;

%Vykreslen� grafick�ho �asov�ho pr�b�hu posloupnosti obdeln�kov�ch impuls�
subplot(4,2,5);
plot(t,w5);
title('Posloupnost posloupnosti obdeln�kov�ch impulsu');
ylabel('p(t)');
xlabel('t');
grid on;

%Vykreslen� grafick�ho �asov�ho pr�b�hu impulsn�ho sign�lu vytvo�en�ho
%p�irozen�m vzorkov�n�m
subplot(4,2,6);
plot(t,w2_p);
title('PAM modulace - p�irozen� vzorkov�n�');
ylabel('w_2_i_p(t)')
xlabel('t');
grid on;

%Vykreslen� grafick�ho �asov�ho pr�b�hu obdeln�kov�ho impulsu
subplot(4,2,7);
plot(t,w6);
title('Obdeln�kov� impuls');
ylabel('p_1(t)');
xlabel('t');
grid on;

%Vykreslen� grafick�ho �asov�ho pr�b�hu impulsn�ho sign�lu vytvo�en�ho
%okam�it�m vzorkov�n�m
subplot(4,2,8);
plot(t,w2_o);
title('PAM modulace - okam�it� vzorkov�n�');
ylabel('w_2_i_o(t)')
xlabel('t');
grid on;

%--------------------------------------------------------------------------
% GRAFY - ZOBRAZEN� PAM MODULACE KONSTANTN�HO SIGN�LU V �ASOV� OBLASTI
%--------------------------------------------------------------------------
% Vykreslen� grafick�ho �asov�ho pr�b�hu konstantn�ho sign�lu
figure;
subplot(4,1,1)
plot(t,w3); 
title('Konstantn� sign�l w_3(t)');
ylabel('w_3(t)')
xlabel('t');
grid on;

%Vykreslen� grafick�ho �asov�ho pr�b�hu posloupnosti Diracov�ch impulsu
subplot(4,2,3);
plot(t,w4);
title('Posloupnost Diracov�ch impuls�');
ylabel('\delta(t)');
xlabel('t');
grid on;

%Vykreslen� grafick�ho �asov�ho pr�b�hu impulsn�ho sign�lu vytvo�en�ho
%ide�ln�m vzorkov�n�m
subplot(4,2,4);
plot(t,w3_i);
title('PAM modulace - ide�ln� vzorkov�n�');
ylabel('w_3_i(t)')
xlabel('t');
grid on;

%Vykreslen� grafick�ho �asov�ho pr�b�hu posloupnosti obdeln�kov�ch impuls�
subplot(4,2,5);
plot(t,w5);
title('Posloupnost posloupnosti obdeln�kov�ch impulsu');
ylabel('p(t)');
xlabel('t');
grid on;

%Vykreslen� grafick�ho �asov�ho pr�b�hu impulsn�ho sign�lu vytvo�en�ho
%p�irozen�m vzorkov�n�m
subplot(4,2,6);
plot(t,w3_p);
title('PAM modulace - p�irozen� vzorkov�n�');
ylabel('w_3_i_p(t)')
xlabel('t');
grid on;

%Vykreslen� grafick�ho �asov�ho pr�b�hu obdeln�kov�ho impulsu
subplot(4,2,7);
plot(t,w6);
title('Obdeln�kov� impuls');
ylabel('p_1(t)');
xlabel('t');
grid on;

%Vykreslen� grafick�ho �asov�ho pr�b�hu impulsn�ho sign�lu vytvo�en�ho
%okam�it�m vzorkov�n�m
subplot(4,2,8);
plot(t,w3_o);
title('PAM modulace - okam�it� vzorkov�n�');
ylabel('w_3_i_o(t)')
xlabel('t');
grid on;


