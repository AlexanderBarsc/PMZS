%--------------------------------------------------------------------------
%           V�ukov� podpora p�em�tu Modulovan� sign�ly
%--------------------------------------------------------------------------
%                       Program FM MODULACE 
%--------------------------------------------------------------------------
% 
% V�po�et a grafick� vykreslen� modulovan�ho vysokofrekven�n�ho sign�lu 
% pomoc� analogov� �hlov� FM modulace (Frequency Modulation)
%
% Voliteln� modula�n� (informa�n�) sign�l a vysokofrekven�n� nosn�
% 
% V�B - Technick� univerzita Ostrava
% Fakulta elektrotechniky a informatiky
% Katedra m��ic� a ��dic� techniky
% 17.listopadu 15
% Ostrava - Poruba
% 708 33
%
% Vypracoval Zden�k Mach��ek 2010
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Vy�i�t�n� v�ech prom�nn�ch pam�ti a uzav�en� v�ech grafick�ch oken
clear all;
close all;

%--------------------------------------------------------------------------
% Specifikace �asov�ho pr�b�hu-�as t[s] a nosn� frekvence fc[Hz]RF sign�lu
f_c=100; % nosn� frekvence vysokofrekven�n�ho RF sign�lu
cas = 0.4; % �as pr�b�hu simulace 
pocet_bodu=10; %po�et bod� v �ase b�hem jedn� periody fc
deltat=1/(f_c*pocet_bodu); % zanedbateln� �asov� okam�ik
t=0:deltat:cas-deltat; % �asov� pr�b�h

%--------------------------------------------------------------------------
% Definov�n� a v�po�et modula�n�ho sign�lu m(t)
% m_t je voliteln� harmonick� nebo obdeln�kov� sign�l

%{
%1.Hamonick� sign�l m(t)
Wmax = 4;  % amplituda modula�n�ho sign�lu m(t)
f_m = 2;    % frekvence modula�n�ho sign�lu m(t)
m_t=Wmax*cos(2*pi*f_m*t); % informa�n� (modula�n�) sign�l m(t)
%}

%%{
%2.Obdeln�kov� sign�l m(t)
Wmax = 1.5;  % velikost modula�n�ho obdeln�kov�ho sign�lu m(t)
t_p1 =length(t)-300; % po�et prvk� v �ase pro hodnotu m1(t) = 0
m1_t(1:t_p1)=0;      % prvn� ��st sign�lu m1(t)
t_p2 =length(t)-200; % po�et prvk� v �ase pro hodnotu m2(t) = Wmax
m2_t(1:t_p2)=Wmax;   % druh� ��st sign�lu m2(t)
t_p3 =((length(t)-300)); % po�et prvk� v �ase pro hodnotu m3(t) = 0
m3_t(1:t_p3)=0;      % t�et� ��st sign�lu m3(t)
m_t=[m1_t m2_t m3_t]; % v�sledn� slo�en� modula�n� sign�l m(t)
%}

%--------------------------------------------------------------------------
% FM modulace - definov�n� a v�po�et p�smov�ho sign�lu v(t)
% amplitudova modulacni slozka R_t = 1
% fazova modulacni slozka Theta_t=Df.integr�l{m(t)}

A_c = 2;     % zes�len� amplitudy modulovan�ho sign�lu v(t)
R_t = 1;     % amplitudova modulacni slozka 
%Df = 50;     % 1.Index f�zov� modulace pro harmonick� sign�l m(t)
Df = 200;   % 2.Index f�zov� modulace pro obdeln�kov� sign�l m(t)
omega_c = 2*pi*f_c; % uhlov� kmito�et nosn�ho sign�lu

%1.V�po�et integr�lu pro harmonick� sign�l m(t)
%Theta_t = Df*(Wmax/(2*pi*f_m)*sin((2*pi*f_m*t))); % F�zov� modula�n� slo�ka

%%{
%2.V�po�et integr�lu pro obdeln�kov� sign�l m(t)
Theta_t1(1:t_p1) = 0;     % prvn� ��st f�zov� modula�n� slo�ky
t2 = 0:deltat:(t_p2*deltat)-deltat;
Theta_t2(1:t_p2) = Df*t2; % druh� ��st f�zov� modula�n� slo�ky
Theta_t3(1:t_p3) = Df*t_p2;     % t�et� ��st f�zov� modula�n� slo�ky
Theta_t=[Theta_t1 Theta_t2 Theta_t3]; % v�sledn� slo�en� f�zov� modula�n� slo�ka Theta(t)
%}

v_t = A_c*R_t*cos(omega_c*t+Theta_t); % p�smov� (modulovan�) sign�l v(t)

%--------------------------------------------------------------------------
%V�po�et frekven�n�ho spektra p�smov�ho sign�lu pomoc� funkce algoritmu FFT
N=length(v_t); % po�et hodnot p�smov�ho sign�lu v(t) 
v_f = (fftshift(fft(v_t)))./N; % komplexn� vektor frekven�n�ho spektra v(t) 
k=-N/2:N/2-1; % pomocn� v�po�et symetrick�ho pole osy x
f = k.*pocet_bodu*f_c./N; % v�po�et x-ov� osy - frekvence f[Hz]

V_amp = abs(v_f);           %amplitudov� frekven�n� spektrum v(t)
V_vykon = V_amp.^2;         %v�konov� frekven�n� spektrum v(t)

% Je-li nulov� hodnota amplitudy na dan� frekvenci je tak� f�ze = 0
for q=1:N
    if ((abs(real(v_f(q))) < 1e-2) && (abs(imag(v_f(q))) < 1e-2)) 
        v_f(q)=0;   
    end;
end;
V_faze = angle(v_f);        %f�zov� frekven�n� spektrum v(t)

%--------------------------------------------------------------------------
% GRAFY - ZOBRAZEN� SIGN�L� V �ASOV� OBLASTI
%--------------------------------------------------------------------------
% Vykreslen� grafick�ho �asov�ho pr�b�hu modula�n�ho sign�lu
figure;
subplot(2,1,1);
plot(t,m_t);
title('Modulacn� sign�l m(t)');
ylabel('m(t)');
xlabel('t[s]');
 
%--------------------------------------------------------------------------
% Vykreslen� grafick�ho �asov�ho pr�b�hu p�smov�ho (modulovan�ho) sign�lu
subplot(2,1,2);
plot(t,v_t);
title('FM modulace - P�smov� sign�l v(t)');
ylabel('v(t)');
xlabel('t[s]');

%--------------------------------------------------------------------------
% GRAFY - ZOBRAZEN� SIGN�L� VE FREKVEN�N� OBLASTI
%--------------------------------------------------------------------------
%vykresleni amplitudov�ho frekven�n�ho spektra 
figure;
subplot(3,1,1);
stem(f,V_amp)
xlabel('f[Hz]')
ylabel('^F^R|W_m|')
grid on;
title('FM modulace - Amplitudov� frekvencni spektrum pasmoveho (modulovaneho)signalu v(t)')
%--------------------------------------------------------------------------
%vykresleni f�zov�ho frekven�n�ho spektra
subplot(3,1,2);
stem(f,V_faze)
xlabel('f[Hz]')
ylabel('\Theta_m')
grid on
title('FM modulace - Fazove frekvencni spektrum pasmoveho (modulovaneho)signalu v(t)')
%--------------------------------------------------------------------------
%vykresleni v�konov�ho frekven�n�ho spektra
subplot(3,1,3);
stem(f,V_vykon)
xlabel('f[Hz]')
ylabel('^F^R|P_m|') 
grid on
title('FM modulace - Vykonove frekvencni spektrum pasmoveho (modulovaneho)signalu v(t)')
