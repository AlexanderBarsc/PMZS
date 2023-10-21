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
% ZADÁNÍ VSTUPNÍHO ANALOGOVÉHO SIGNÁLU - harmonický signal w(t)
%--------------------------------------------------------------------------
W_max = 4; % Amplituda signálu
f = 1/6; % frekvence signálu
omega = 2*pi*f; % úhlový kmitoèet 
fi = pi; % fázový posun
C = 10; % konstantní složka
w = W_max*cos(omega*t+fi) + C; % harmonický signal w(t) 

%--------------------------------------------------------------------------
% VZORKOVÁNÍ SIGNÁLU - OKAMŽITÉ VZORKOVÁNÍ S ŠÍØKOU IMPULSU TÉMÌØ ROVNU PERIODÌ VZORKOVÁNÍ
%--------------------------------------------------------------------------
% Výpoèet posloupnosti Diracových impulsu
t1=0;  
w1=0;

for n=0:T_s:cas-T_s+2*deltat;
t1a=n:deltat:n;
w1a(1:length(t1a))=1;

t1b=n+deltat:deltat:n+T_s-deltat;
w1b(1:length(t1b))=0;

w1=[w1 w1a w1b]; % signal delta(t)
end;

w_i=w.*w1; %PAM - ideální vzorkování vstupního signálu
sirka = 0.4;  % šíøka impulsu pozor! do max. velikosti menší než T_s
%Definování jednoho impulsu dané šíøky
t2a=0:deltat:sirka;
w2a(1:length(t2a))=1;

t2b=sirka+deltat:deltat:cas;
w2b(1:length(t2b))=0;

w2=[w2a w2b]; % impulsní obdelníkový signal p1(t)

w_o1=conv(w_i,w2); %PAM - okamžité vzorkování vstupního signálu

w_o = w_o1(1:length(t));

%--------------------------------------------------------------------------
% KVANTOVÁNÍ SIGNÁLU - ROVNOMÌRNÉ
%--------------------------------------------------------------------------
amp_rozsah = 15; % pozor! musí být shodná nebo vìtší než hodnota vstupního signálu
pocet_stavu = 4; % pozor! nefunguje v pøípadì vstupního signálu procházejícícho nulou
bit_cislo = 2 % pozor! upravit dle poctu stavu
i1 = 0;
i2 = 1;
%w_kvant1 = zeros(1,length(t));
w_kvant = zeros(1,length(t));
for i = 0:deltat:cas-deltat
    
    if (mod(i,T_s) == 0) || (i == 0)
        i1 = i1+1;
        w_kvant1(i1) = fix(w_o(i2+1)/(amp_rozsah/pocet_stavu));
    end;
    w_kvant(i2) = w_kvant1(i1)*(amp_rozsah/pocet_stavu);
    i2 = i2+1;
end;

%--------------------------------------------------------------------------
% KODOVÁNÍ SIGNÁLU - SERIOVÝ TOK DAT
%--------------------------------------------------------------------------

w_kod1=dec2bin(w_kvant1-1);
n2 = 0;
n3 = 1;
n4 = 1;
for n1 = 0:deltat:cas-deltat
    if (mod(n1,T_s/bit_cislo) == 0) || (n1 == 0)
        n2 = n2+1;
        if n2 > bit_cislo
            n3 = n3+1;
            n2 = 1;
        end;
    end;
    w_kod(n4)=w_kod1(n3,n2);
    n4 = n4+1;
end;

w_kodovany = (double(w_kod))-48;



%--------------------------------------------------------------------------
% GRAFY - ZOBRAZENÍ PCM MODULACE VSTUPNÍHO SIGNÁLU V ÈASOVÉ OBLASTI
%--------------------------------------------------------------------------
% Vykreslení grafického èasového prùbìhu vstupního signálu
figure;
subplot(4,1,1)
plot(t,w); 
title('Vstupní signál w(t)');
ylabel('w(t)')
xlabel('t');
grid on;

%Vykreslení grafického èasového prùbìhu impulsního signálu vytvoøeného
%okamžitým vzorkováním
subplot(4,1,2);
plot(t,w_o);
title('VZORKOVÁNÍ - PAM modulace -  okamžité vzorkování');
ylabel('w_i_o(t)')
xlabel('t');
grid on;

%Vykreslení grafického èasového prùbìhu impulsního signálu vytvoøeného
%okamžitým vzorkováním
subplot(4,1,3);
plot(t,w_kvant);
title('KVANTOVÁNÍ - rovnomìrné kvantování');
ylabel('w_k_v_a_n_t(t)')
xlabel('t');
grid on;
%Vykreslení grafického èasového prùbìhu impulsního signálu vytvoøeného
%okamžitým vzorkováním
subplot(4,1,4);
plot(t(1:(length(t)-1)),w_kodovany);
title('KÓDOVÁNÍ - kodovaný tok dat z kvantovaného signálu - BEZ dalsího kodovaní');
ylabel('w_k_o_d_o_v_a_n_y(t)')
xlabel('t');
grid on;


