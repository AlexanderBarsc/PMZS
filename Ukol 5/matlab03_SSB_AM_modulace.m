%--------------------------------------------------------------------------
%           V�ukov� podpora p�em�tu Modulovan� sign�ly
%--------------------------------------------------------------------------
%                       Program SSB-AM MODULACE (USSB-AM, LSSB-AM)
%--------------------------------------------------------------------------
% 
% V�po�et a grafick� vykreslen� modulovan�ho vysokofrekven�n�ho sign�lu 
% pomoc� analogov� amplitudov� SSB-AM
% (SingleSide Band - Amplitude Modulation)
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
f_c=1000; % nosn� frekvence vysokofrekven�n�ho RF sign�lu
cas = 0.1; % �as pr�b�hu simulace 
pocet_bodu=3; %po�et bod� v �ase b�hem jedn� periody fc
deltat=1/(f_c*pocet_bodu); % zanedbateln� �asov� okam�ik
t=0:deltat:cas-deltat; % �asov� pr�b�h

%--------------------------------------------------------------------------
% Definov�n� a v�po�et modula�n�ho sign�lu m(t)
% m_t je slo�eny ze dvou harmonick�ch slo�ek

Wmax1 = 3;    % amplituda 1.harmonick� slo�ky modula�n�ho sign�lu m(t)
f_m1 = 100;   % frekvence 1.harmonick� slo�ky modula�n�ho sign�lu m(t)
faze1 = pi/8; % f�ze 1.harmonick� slo�ky modula�n�ho sign�lu m(t)
Wmax2 = 5;    % amplituda 2.harmonick� slo�ky modula�n�ho sign�lu m(t)
f_m2 = 40;    % frekvence 2.harmonick� slo�ky modula�n�ho sign�lu m(t)
faze2 = pi;   % f�ze 2.harmonick� slo�ky modula�n�ho sign�lu m(t)
% informa�n� (modula�n�) sign�l m(t)
m_t=Wmax1*cos(2*pi*f_m1*t+faze1)+Wmax2*cos(2*pi*f_m2*t+faze2); 
%--------------------------------------------------------------------------
% v�po�et Hilbertovy transformace signalu m(t)
m_H_t=Wmax1*sin(2*pi*f_m1*t+faze1)+Wmax2*sin(2*pi*f_m2*t+faze2); 

%--------------------------------------------------------------------------
% SSB-AM modulace - definov�n� a v�po�et p�smov�ho sign�lu v(t)
% komplexn� ob�lka g(t) = m(t)+-jmH(t)
% v(t)= Re(g(t)).exp(j.omega.fc.t)

% 1.USSB-AM modulace - Upper Single Side Band - Amplitude Modulation
%v_t = real((m_t+j*m_H_t).*exp(j*2*pi*f_c*t)); % p�smov� (modulovan�) sign�l v(t)

% 2.LSSB-AM modulace - Lower Single Side Band - Amplitude Modulation
v_t = real((m_t-j*m_H_t).*exp(j*2*pi*f_c*t)); % p�smov� (modulovan�) sign�l v(t)

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
    if ((abs(real(v_f(q))) < 1e-3) && (abs(imag(v_f(q))) < 1e-3)) 
        v_f(q)=0;   
    end;
end;
V_faze = angle(v_f);        %f�zov� frekven�n� spektrum v(t)

%--------------------------------------------------------------------------
% GRAFY - ZOBRAZEN� SIGN�L� V �ASOV� OBLASTI
%--------------------------------------------------------------------------
% Vykreslen� grafick�ho �asov�ho pr�b�hu modula�n�ho sign�lu
figure;
subplot(3,1,1);
plot(t,m_t);
title('Modulacn� sign�l m(t)');
ylabel('m(t)');
xlabel('t[s]');
 
%--------------------------------------------------------------------------
% Vykreslen� grafick�ho �asov�ho pr�b�hu transformovan�ho modula�n�ho sign�lu
% po Hilbertov� transformaci
subplot(3,1,2);
plot(t,m_H_t);
title('Transformovan� modulacn� sign�l po Hilbertove transformaci m_H(t)');
ylabel('m_H(t)');
xlabel('t[s]');

%--------------------------------------------------------------------------
% Vykreslen� grafick�ho �asov�ho pr�b�hu p�smov�ho (modulovan�ho) sign�lu
subplot(3,1,3);
plot(t,v_t);
title('SSB-AM modulace - P�smov� sign�l v(t)');
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
title('SSB-AM modulace - Amplitudov� frekvencni spektrum pasmoveho (modulovaneho)signalu v(t)')
%--------------------------------------------------------------------------
%vykresleni f�zov�ho frekven�n�ho spektra
subplot(3,1,2);
stem(f,V_faze)
xlabel('f[Hz]')
ylabel('\Theta_m')
grid on
title('SSB-AM modulace - Fazove frekvencni spektrum pasmoveho (modulovaneho)signalu v(t)')
%--------------------------------------------------------------------------
%vykresleni v�konov�ho frekven�n�ho spektra
subplot(3,1,3);
stem(f,V_vykon)
xlabel('f[Hz]')
ylabel('^F^R|P_m|') 
grid on
title('SSB-AM modulace - Vykonove frekvencni spektrum pasmoveho (modulovaneho)signalu v(t)')

