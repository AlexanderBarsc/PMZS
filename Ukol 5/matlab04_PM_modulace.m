%--------------------------------------------------------------------------
%           V�ukov� podpora p�em�tu Modulovan� sign�ly
%--------------------------------------------------------------------------
%                       Program PM MODULACE 
%--------------------------------------------------------------------------
% 
% V�po�et a grafick� vykreslen� modulovan�ho vysokofrekven�n�ho sign�lu 
% pomoc� analogov� �hlov� PM modulace (Phase Modulation)
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
f_c=200; % nosn� frekvence vysokofrekven�n�ho RF sign�lu
cas = 0.2; % �as pr�b�hu simulace 
pocet_bodu=10; %po�et bod� v �ase b�hem jedn� periody fc
deltat=1/(f_c*pocet_bodu); % zanedbateln� �asov� okam�ik
t=0:deltat:cas-deltat; % �asov� pr�b�h

%--------------------------------------------------------------------------
% Definov�n� a v�po�et modula�n�ho sign�lu m(t)
% m_t je voliteln� harmonick� nebo obdeln�kov� sign�l

%{
%1.Hamonick� sign�l m(t)
Wmax = 1.5;  % amplituda modula�n�ho sign�lu m(t)
f_m = 10;    % frekvence modula�n�ho sign�lu m(t)
faze = pi/2; % f�ze modula�n�ho sign�lu m(t)
m_t=Wmax*cos(2*pi*f_m*t+faze); % informa�n� (modula�n�) sign�l m(t)
%}

%2.Obdeln�kov� sign�l m(t)
Wmax = 1.5;  % velikost modula�n�ho obdeln�kov�ho sign�lu m(t)
t_p1 =length(t)-300; % po�et prvk� v �ase pro hodnotu m1(t) = 0
m1_t(1:t_p1)=0;      % prvn� ��st sign�lu m1(t)
t_p2 =length(t)-200; % po�et prvk� v �ase pro hodnotu m2(t) = Wmax
m2_t(1:t_p2)=Wmax;   % druh� ��st sign�lu m2(t)
t_p3 =((length(t)-300)); % po�et prvk� v �ase pro hodnotu m3(t) = 0
m3_t(1:t_p3)=0;      % t�et� ��st sign�lu m3(t)
m_t=[m1_t m2_t m3_t]; % v�sledn� slo�en� modula�n� sign�l m(t)


%--------------------------------------------------------------------------
% PM modulace - definov�n� a v�po�et p�smov�ho sign�lu v(t)
% amplitudova modulacni slozka R_t = 1
% fazova modulacni slozka Theta_t=Dp.m(t)

A_c = 2;            % zes�len� amplitudy modulovan�ho sign�lu v(t)
R_t=1;              % amplitudova modulacni slozka 
Dp = 1;             % Index f�zov� modulace
Theta_t = Dp*m_t;   % F�zov� modula�n� slo�ka
omega_c = 2*pi*f_c; % uhlov� kmito�et nosn�ho sign�lu
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
    if ((abs(real(v_f(q))) < 3e-2) && (abs(imag(v_f(q))) < 3e-2)) 
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
title('PM modulace - P�smov� sign�l v(t)');
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
title('PM modulace - Amplitudov� frekvencni spektrum pasmoveho (modulovaneho)signalu v(t)')
%--------------------------------------------------------------------------
%vykresleni f�zov�ho frekven�n�ho spektra
subplot(3,1,2);
stem(f,V_faze)
xlabel('f[Hz]')
ylabel('\Theta_m')
grid on
title('PM modulace - Fazove frekvencni spektrum pasmoveho (modulovaneho)signalu v(t)')
%--------------------------------------------------------------------------
%vykresleni v�konov�ho frekven�n�ho spektra
subplot(3,1,3);
stem(f,V_vykon)
xlabel('f[Hz]')
ylabel('^F^R|P_m|') 
grid on
title('PM modulace - Vykonove frekvencni spektrum pasmoveho (modulovaneho)signalu v(t)')
