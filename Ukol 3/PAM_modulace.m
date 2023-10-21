%--------------------------------------------------------------------------
% Vyèištìní všech promìnných pamìti a uzavøení všech grafických oken
clear all;
close all;

%--------------------------------------------------------------------------
% Specifikace èasového prùbìhu-èas t[s] signálu
cas = 10; % èas prùbìhu simulace 
pocet_bodu = 100; %poèet bodù v èase bìhem jedné vzorkovací periody fs
T_s = 0.5; % Perioda vzorkování
deltat = T_s/pocet_bodu; % zanedbatelný èasový okamžik
t = 0:deltat:cas; % èasový prùbìh
%--------------------------------------------------------------------------
% 1. zadání harmonického signalu w1(t)
W_max = 10; % Amplituda signálu
f = 1/4; % frekvence signálu
omega = 2*pi*f; % úhlový kmitoèet 
fi = pi; % fázový posun
C = 1; % konstantní složka
w1 = W_max*cos(omega*t+fi) + C; % harmonický signal w1(t) 
%--------------------------------------------------------------------------
% 2. zadání obdelníkového signalu w2(t)
t2a=0:deltat:2-deltat;
w2a(1:length(t2a))=0; %1.èást signalu w2(t)

t2b=2:deltat:4;
w2b(1:length(t2b))=3; %2.èást signalu w2(t)

t2c=4+deltat:deltat:cas;
w2c(1:length(t2c))=0; %3.èást signalu w2(t)

w2=[w2a w2b w2c]; % složený celkový signal w2(t)
%--------------------------------------------------------------------------
% 3. zadání konstantního signálu w3(t)
w3(1:length(t))=5; %konstantní signal w3(t) 

%--------------------------------------------------------------------------
% VÝPOÈET IMPULSNÍHO SIGNÁLU - IDEÁLNÍ VZORKOVÁNÍ
%--------------------------------------------------------------------------
% Výpoèet posloupnosti Diracových impulsu
t4=0;  
w4=0;

for n=0:T_s:cas-T_s+2*deltat;
t4a=n:deltat:n;
w4a(1:length(t4a))=1;

t4b=n+deltat:deltat:n+T_s-deltat;
w4b(1:length(t4b))=0;

w4=[w4 w4a w4b]; % signal delta(t)
end;

w1_i=w1.*w4; %PAM - ideální vzorkování harmonického signálu
w2_i=w2.*w4; %PAM - ideální vzorkování obdelníkového signálu
w3_i=w3.*w4; %PAM - ideální vzorkování konstantního signálu

%--------------------------------------------------------------------------
% VÝPOÈET IMPULSNÍHO SIGNÁLU - PØIROZENÉ VZORKOVÁNÍ
%--------------------------------------------------------------------------
%Definování umpulsù dané šíøky
sirka = 0.2;  % šíøka impulsu
t5=0;
w5=0;

for n=0:T_s:cas-T_s+2*deltat;
t5a=n:deltat:n+sirka;
w5a(1:length(t5a))=1;

t5b=n+sirka+deltat:deltat:n+T_s-deltat;
w5b(1:length(t5b))=0;

w5=[w5 w5a w5b]; % posloupnost impulsních obdelníkových signalu p(t)
end;

w1_p=w1.*w5; %PAM - pøirozené vzorkování harmonického signálu
w2_p=w2.*w5; %PAM - pøirozené vzorkování obdelníkového signálu
w3_p=w3.*w5; %PAM - pøirozené vzorkování konstantního signálu

%--------------------------------------------------------------------------
% VÝPOÈET IMPULSNÍHO SIGNÁLU - OKAMŽITÉ VZORKOVÁNÍ
%--------------------------------------------------------------------------
%Definování jednoho impulsu dané šíøky
t6a=0:deltat:sirka;
w6a(1:length(t6a))=1;

t6b=sirka+deltat:deltat:cas;
w6b(1:length(t6b))=0;

w6=[w6a w6b]; % impulsní obdelníkový signal p1(t)

w1_o1=conv(w1_i,w6); %PAM - okamžité vzorkování harmonického signálu
w1_o = w1_o1(1:length(t));
w2_o1=conv(w2_i,w6); %PAM - okamžité vzorkování obdelníkového signálu
w2_o = w2_o1(1:length(t));
w3_o1=conv(w3_i,w6); %PAM - okamžité vzorkování konstantního signálu
w3_o = w3_o1(1:length(t));

%--------------------------------------------------------------------------
% GRAFY - ZOBRAZENÍ PAM MODULACE HARMONICKÉHO SIGNÁLU V ÈASOVÉ OBLASTI
%--------------------------------------------------------------------------
% Vykreslení grafického èasového prùbìhu harmonického signálu
figure;
subplot(4,1,1)
plot(t,w1); 
title('Harmonický signál w_1(t)');
ylabel('w_1(t)')
xlabel('t');
grid on;

%Vykreslení grafického èasového prùbìhu posloupnosti Diracových impulsu
subplot(4,2,3);
plot(t,w4);
title('Posloupnost Diracových impulsù');
ylabel('\delta(t)');
xlabel('t');
grid on;

%Vykreslení grafického èasového prùbìhu impulsního signálu vytvoøeného
%ideálním vzorkováním
subplot(4,2,4);
plot(t,w1_i);
title('PAM modulace - ideální vzorkování');
ylabel('w_1_i(t)')
xlabel('t');
grid on;

%Vykreslení grafického èasového prùbìhu posloupnosti obdelníkových impulsù
subplot(4,2,5);
plot(t,w5);
title('Posloupnost posloupnosti obdelníkových impulsu');
ylabel('p(t)');
xlabel('t');
grid on;

%Vykreslení grafického èasového prùbìhu impulsního signálu vytvoøeného
%pøirozeným vzorkováním
subplot(4,2,6);
plot(t,w1_p);
title('PAM modulace - pøirozené vzorkování');
ylabel('w_1_i_p(t)')
xlabel('t');
grid on;

%Vykreslení grafického èasového prùbìhu obdelníkového impulsu
subplot(4,2,7);
plot(t,w6);
title('Obdelníkový impuls');
ylabel('p_1(t)');
xlabel('t');
grid on;

%Vykreslení grafického èasového prùbìhu impulsního signálu vytvoøeného
%okamžitým vzorkováním
subplot(4,2,8);
plot(t,w1_o);
title('PAM modulace - okamžité vzorkování');
ylabel('w_1_i_o(t)')
xlabel('t');
grid on;

%--------------------------------------------------------------------------
% GRAFY - ZOBRAZENÍ PAM MODULACE PERIODICKÉHO OBDELNÍKOVÉHO SIGNÁLU V ÈASOVÉ OBLASTI
%--------------------------------------------------------------------------
% Vykreslení grafického èasového prùbìhu periodického obdelníkového signálu
figure;
subplot(4,1,1)
plot(t,w2); 
title('Obdelníkový periodický signál w_2(t)');
ylabel('w_2(t)')
xlabel('t');
grid on;

%Vykreslení grafického èasového prùbìhu posloupnosti Diracových impulsu
subplot(4,2,3);
plot(t,w4);
title('Posloupnost Diracových impulsù');
ylabel('\delta(t)');
xlabel('t');
grid on;

%Vykreslení grafického èasového prùbìhu impulsního signálu vytvoøeného
%ideálním vzorkováním
subplot(4,2,4);
plot(t,w2_i);
title('PAM modulace - ideální vzorkování');
ylabel('w_2_i(t)')
xlabel('t');
grid on;

%Vykreslení grafického èasového prùbìhu posloupnosti obdelníkových impulsù
subplot(4,2,5);
plot(t,w5);
title('Posloupnost posloupnosti obdelníkových impulsu');
ylabel('p(t)');
xlabel('t');
grid on;

%Vykreslení grafického èasového prùbìhu impulsního signálu vytvoøeného
%pøirozeným vzorkováním
subplot(4,2,6);
plot(t,w2_p);
title('PAM modulace - pøirozené vzorkování');
ylabel('w_2_i_p(t)')
xlabel('t');
grid on;

%Vykreslení grafického èasového prùbìhu obdelníkového impulsu
subplot(4,2,7);
plot(t,w6);
title('Obdelníkový impuls');
ylabel('p_1(t)');
xlabel('t');
grid on;

%Vykreslení grafického èasového prùbìhu impulsního signálu vytvoøeného
%okamžitým vzorkováním
subplot(4,2,8);
plot(t,w2_o);
title('PAM modulace - okamžité vzorkování');
ylabel('w_2_i_o(t)')
xlabel('t');
grid on;

%--------------------------------------------------------------------------
% GRAFY - ZOBRAZENÍ PAM MODULACE KONSTANTNÍHO SIGNÁLU V ÈASOVÉ OBLASTI
%--------------------------------------------------------------------------
% Vykreslení grafického èasového prùbìhu konstantního signálu
figure;
subplot(4,1,1)
plot(t,w3); 
title('Konstantní signál w_3(t)');
ylabel('w_3(t)')
xlabel('t');
grid on;

%Vykreslení grafického èasového prùbìhu posloupnosti Diracových impulsu
subplot(4,2,3);
plot(t,w4);
title('Posloupnost Diracových impulsù');
ylabel('\delta(t)');
xlabel('t');
grid on;

%Vykreslení grafického èasového prùbìhu impulsního signálu vytvoøeného
%ideálním vzorkováním
subplot(4,2,4);
plot(t,w3_i);
title('PAM modulace - ideální vzorkování');
ylabel('w_3_i(t)')
xlabel('t');
grid on;

%Vykreslení grafického èasového prùbìhu posloupnosti obdelníkových impulsù
subplot(4,2,5);
plot(t,w5);
title('Posloupnost posloupnosti obdelníkových impulsu');
ylabel('p(t)');
xlabel('t');
grid on;

%Vykreslení grafického èasového prùbìhu impulsního signálu vytvoøeného
%pøirozeným vzorkováním
subplot(4,2,6);
plot(t,w3_p);
title('PAM modulace - pøirozené vzorkování');
ylabel('w_3_i_p(t)')
xlabel('t');
grid on;

%Vykreslení grafického èasového prùbìhu obdelníkového impulsu
subplot(4,2,7);
plot(t,w6);
title('Obdelníkový impuls');
ylabel('p_1(t)');
xlabel('t');
grid on;

%Vykreslení grafického èasového prùbìhu impulsního signálu vytvoøeného
%okamžitým vzorkováním
subplot(4,2,8);
plot(t,w3_o);
title('PAM modulace - okamžité vzorkování');
ylabel('w_3_i_o(t)')
xlabel('t');
grid on;


