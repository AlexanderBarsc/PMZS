%--------------------------------------------------------------------------
%           V�ukov� podpora p�em�tu Modulovan� sign�ly
%--------------------------------------------------------------------------
%                       Program QM MODULACE 
%--------------------------------------------------------------------------
% 
% V�po�et a grafick� vykreslen� modulovan�ho vysokofrekven�n�ho sign�lu 
% pomoc� analogov� kvadraturn� QM modulace (Quadrature Modulation)
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
f_c=5000; % nosn� frekvence vysokofrekven�n�ho RF sign�lu
cas = 0.05; % �as pr�b�hu simulace 
pocet_bodu=10; %po�et bod� v �ase b�hem jedn� periody fc
deltat=1/(f_c*pocet_bodu); % zanedbateln� �asov� okam�ik
t=0:deltat:cas-deltat; % �asov� pr�b�h

%--------------------------------------------------------------------------
% Definov�n� a v�po�et modula�n�ho sign�lu m(t)
% m_t je voliteln� harmonick� nebo obdeln�kov� sign�l

%1.Hamonick� sign�l m(t)
Wmax1 = 4;  % amplituda modula�n�ho sign�lu m(t)
f_m1 = 100;    % frekvence modula�n�ho sign�lu m(t)
m_t1 = Wmax1*cos(2*pi*f_m1*t); % informa�n� (modula�n�) sign�l m(t)

%2.Hamonick� sign�l m(t)
Wmax2 = 6;  % amplituda modula�n�ho sign�lu m(t)
f_m2 = 400;    % frekvence modula�n�ho sign�lu m(t)
m_t2 = Wmax2*cos(2*pi*f_m2*t); % informa�n� (modula�n�) sign�l m(t)


%--------------------------------------------------------------------------
% QM modulace - definov�n� a v�po�et p�smov�ho sign�lu v(t)
omega_c = 2*pi*f_c; % uhlov� kmito�et nosn�ho sign�lu

%V�po�et komplexn� ob�lky sign�lu v(t)=m1(t)+j.m2(t)
g_t = m_t1 + 1j*m_t2;          % komplexni obalka

%Vypocet pasmoveho (modulovaneho) signalu z algoritmu komplexn� ob�lky
v_t = real(g_t.*exp(1j*omega_c*t));

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
% Vykreslen� grafick�ho �asov�ho pr�b�hu modula�n�ho sign�lu m1(t)
figure;
subplot(3,1,1);
plot(t,m_t1);
title('Modulacn� sign�l m_1(t)');
ylabel('m_1(t)');
xlabel('t[s]');
 
%--------------------------------------------------------------------------
% Vykreslen� grafick�ho �asov�ho pr�b�hu modula�n�ho sign�lu m2(t)
subplot(3,1,2);
plot(t,m_t2);
title('Modulacn� sign�l m_2(t)');
ylabel('m_2(t)');
xlabel('t[s]');

%--------------------------------------------------------------------------
% Vykreslen� grafick�ho �asov�ho pr�b�hu p�smov�ho (modulovan�ho) sign�lu
subplot(3,1,3);
plot(t,v_t);
title('QM modulace - P�smov� sign�l v(t)');
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
title('QM modulace - Amplitudov� frekvencni spektrum pasmoveho (modulovaneho)signalu v(t)')
%--------------------------------------------------------------------------
%vykresleni f�zov�ho frekven�n�ho spektra
subplot(3,1,2);
stem(f,V_faze)
xlabel('f[Hz]')
ylabel('\Theta_m')
grid on
title('QM modulace - Fazove frekvencni spektrum pasmoveho (modulovaneho)signalu v(t)')
%--------------------------------------------------------------------------
%vykresleni v�konov�ho frekven�n�ho spektra
subplot(3,1,3);
stem(f,V_vykon)
xlabel('f[Hz]')
ylabel('^F^R|P_m|') 
grid on
title('QM modulace - Vykonove frekvencni spektrum pasmoveho (modulovaneho)signalu v(t)')



