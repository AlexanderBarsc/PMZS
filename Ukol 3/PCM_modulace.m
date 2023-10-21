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
% ZAD�N� VSTUPN�HO ANALOGOV�HO SIGN�LU - harmonick� signal w(t)
%--------------------------------------------------------------------------
W_max = 4; % Amplituda sign�lu
f = 1/6; % frekvence sign�lu
omega = 2*pi*f; % �hlov� kmito�et 
fi = pi; % f�zov� posun
C = 10; % konstantn� slo�ka
w = W_max*cos(omega*t+fi) + C; % harmonick� signal w(t) 

%--------------------------------------------------------------------------
% VZORKOV�N� SIGN�LU - OKAM�IT� VZORKOV�N� S ���KOU IMPULSU T�M�� ROVNU PERIOD� VZORKOV�N�
%--------------------------------------------------------------------------
% V�po�et posloupnosti Diracov�ch impulsu
t1=0;  
w1=0;

for n=0:T_s:cas-T_s+2*deltat;
t1a=n:deltat:n;
w1a(1:length(t1a))=1;

t1b=n+deltat:deltat:n+T_s-deltat;
w1b(1:length(t1b))=0;

w1=[w1 w1a w1b]; % signal delta(t)
end;

w_i=w.*w1; %PAM - ide�ln� vzorkov�n� vstupn�ho sign�lu
sirka = 0.4;  % ���ka impulsu pozor! do max. velikosti men�� ne� T_s
%Definov�n� jednoho impulsu dan� ���ky
t2a=0:deltat:sirka;
w2a(1:length(t2a))=1;

t2b=sirka+deltat:deltat:cas;
w2b(1:length(t2b))=0;

w2=[w2a w2b]; % impulsn� obdeln�kov� signal p1(t)

w_o1=conv(w_i,w2); %PAM - okam�it� vzorkov�n� vstupn�ho sign�lu

w_o = w_o1(1:length(t));

%--------------------------------------------------------------------------
% KVANTOV�N� SIGN�LU - ROVNOM�RN�
%--------------------------------------------------------------------------
amp_rozsah = 15; % pozor! mus� b�t shodn� nebo v�t�� ne� hodnota vstupn�ho sign�lu
pocet_stavu = 4; % pozor! nefunguje v p��pad� vstupn�ho sign�lu proch�zej�c�cho nulou
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
% KODOV�N� SIGN�LU - SERIOV� TOK DAT
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
% GRAFY - ZOBRAZEN� PCM MODULACE VSTUPN�HO SIGN�LU V �ASOV� OBLASTI
%--------------------------------------------------------------------------
% Vykreslen� grafick�ho �asov�ho pr�b�hu vstupn�ho sign�lu
figure;
subplot(4,1,1)
plot(t,w); 
title('Vstupn� sign�l w(t)');
ylabel('w(t)')
xlabel('t');
grid on;

%Vykreslen� grafick�ho �asov�ho pr�b�hu impulsn�ho sign�lu vytvo�en�ho
%okam�it�m vzorkov�n�m
subplot(4,1,2);
plot(t,w_o);
title('VZORKOV�N� - PAM modulace -  okam�it� vzorkov�n�');
ylabel('w_i_o(t)')
xlabel('t');
grid on;

%Vykreslen� grafick�ho �asov�ho pr�b�hu impulsn�ho sign�lu vytvo�en�ho
%okam�it�m vzorkov�n�m
subplot(4,1,3);
plot(t,w_kvant);
title('KVANTOV�N� - rovnom�rn� kvantov�n�');
ylabel('w_k_v_a_n_t(t)')
xlabel('t');
grid on;
%Vykreslen� grafick�ho �asov�ho pr�b�hu impulsn�ho sign�lu vytvo�en�ho
%okam�it�m vzorkov�n�m
subplot(4,1,4);
plot(t(1:(length(t)-1)),w_kodovany);
title('K�DOV�N� - kodovan� tok dat z kvantovan�ho sign�lu - BEZ dals�ho kodovan�');
ylabel('w_k_o_d_o_v_a_n_y(t)')
xlabel('t');
grid on;


