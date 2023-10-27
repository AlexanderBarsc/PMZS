%% PMZS Cv.04

clear all
close all

%% Definice zadaneho harmonickeho signalu

amplitude = 0.2;         
omega = 10*pi;
phase = -pi/2;
f = omega/(2*pi);
Tm = 1/f;

% Specifikace časového průběhu-čas t[s] a nosné frekvence fc[Hz]RF signálu
f_c = 50; % nosná frekvence vysokofrekvenčního RF signálu (zvoleno)
cas = 2*Tm; % čas průběhu simulace (2 periody)
pocet_bodu = 3; %počet bodů v čase během jedné periody fc (zvoleno)
deltat = 1/(f_c*pocet_bodu); % zanedbatelný časový okamžik

t = 0:deltat:cas-deltat; % časový průběh

m_t = amplitude*cos(omega*t+phase);		% signal m[t]

% Vykresleni signalu
figure()
plot(t,m_t); 
xlabel('t(s)')
ylabel('m(t)')
title('Harmonický průběh m(t)')

%% Amplitudove a fazove spektrum

% Provedeni FFT
m_tt0=fft(m_t);
m_t0=fftshift(m_tt0);

k = (-length(m_t0)/2) :1: (length(m_t0)/2)-1;
N = length(m_t0);
fm = k.*pocet_bodu*f_c./N;

figure;
subplot(2,1,1)
stem(fm,(abs(m_t0)/length(m_t0)))
xlabel('f(Hz)')
ylabel('|m(f)|')
title('Amplitudové spektrum |m(f)|')

subplot(2,1,2)
stem(fm,(atan2(round(imag(m_t0),6),round(real(m_t0),6))))
xlabel('f(Hz)')
ylabel('phase m(f)')
title('Fázové spektrum m(f)')

%% AM MODULACE

% amplitudova modulacni slozka R_t = |1+m_t|
% fazova modulacni slozka Theta_t=0 pro m(t)>-1, Theta_t=pi pro m(t)<-1;
% m(t)nesmi byt premodulovany >1

R_t=abs(1+m_t); % amplitudova modulacni slozka 
if (m_t>-1)     % fázová modulační složka
    Theta_t = 0;
else 
    Theta_t = pi;
end
v_am_t=R_t.*cos(2*pi*f_c*t+Theta_t); % pásmový (modulovaný) signál v(t)

%Výpočet frekvenčního spektra pásmového signálu pomocí funkce algoritmu FFT
N_am = length(v_am_t);                                    % počet hodnot pásmového signálu v(t) 
v_am_f = (fftshift(fft(v_am_t)))./N_am;                   % komplexní vektor frekvenčního spektra v(t) 
k_am = -N_am/2:N_am/2-1;                                  % pomocný výpočet symetrického pole osy x
f_am = k_am.*pocet_bodu*f_c./N_am;                        % výpočet x-ové osy - frekvence f[Hz]

V_am_amp = abs(v_am_f);                                   % amplitudové frekvenční spektrum v(t)

for q=1:N_am
    if ((abs(real(v_am_f(q))) < 1e-3) && (abs(imag(v_am_f(q))) < 1e-3)) 
        v_am_f(q)=0;                                   % Je-li nulová hodnota amplitudy na dané frekvenci je také fáze = 0
    end;
end;

V_am_faze = angle(v_am_f);                             %fázové frekvenční spektrum v(t)

% Vykreslení grafického časového průběhu pásmového signálu
figure;
plot(t,v_am_t);
title('AM modulace - Pásmový signál v(t)');
ylabel('v(t)');
xlabel('t(s)');

%vykresleni amplitudového frekvenčního spektra 
figure;
subplot(2,1,1);
stem(f_am,V_am_amp)
xlabel('f(Hz)')
ylabel('^F^R|W_m|')
grid on;
title('AM modulace - Amplitudové frekvencni spektrum pasmoveho signalu v(t)')

%vykresleni fázového frekvenčního spektra
subplot(2,1,2);
stem(f_am,V_am_faze)
xlabel('f(Hz)')
ylabel('\Theta_m')
grid on
title('AM modulace - Fazove frekvencni spektrum pasmoveho signalu v(t)')

%% DSB_SC MODULACE

% amplitudova modulacni slozka R_t = |m_t|
% fazova modulacni slozka Theta_t=0 pro m(t)>0, Theta_t=pi pro m(t)<0;

x_t = m_t; 											% soufázová modulacni slozka
y_t = 0;   											% kvadraturní modulacni slozka

v_db_t = x_t.*cos(2*pi*f_c*t)-y_t.*sin(2*pi*f_c*t); % pásmový (modulovaný) signál v(t)

%Výpočet frekvenčního spektra pásmového signálu pomocí funkce algoritmu FFT
N_db=length(v_db_t);                                % počet hodnot pásmového signálu v(t) 
v_db_f = (fftshift(fft(v_db_t)))./N_db;             % komplexní vektor frekvenčního spektra v(t) 
k_db=-N_db/2:N_db/2-1;                              % pomocný výpočet symetrického pole osy x
f_db = k_db.*pocet_bodu*f_c./N_db;                  % výpočet x-ové osy - frekvence f[Hz]

V_db_amp = abs(v_db_f);                             %amplitudové frekvenční spektrum v(t)

for q=1:N_db                                        
    if ((abs(real(v_db_f(q))) < 1e-3) && (abs(imag(v_db_f(q))) < 1e-3))  
        v_db_f(q)=0;                                % Je-li nulová hodnota amplitudy na dané frekvenci je také fáze = 0
    end;
end;

V_db_faze = angle(v_db_f);                           %fázové frekvenční spektrum v(t)

% Vykreslení grafického časového průběhu pásmového (modulovaného) signálu
figure;
plot(t,v_db_t);
title('DSB-SC modulace - Pásmový signál v(t)');
ylabel('v(t)');
xlabel('t(s)');

%vykresleni amplitudového frekvenčního spektra 
figure;
subplot(2,1,1);
stem(f_db,V_db_amp)
xlabel('f(Hz)')
ylabel('^F^R|W_m|')
grid on;
title('DSB-SC modulace - Amplitudové kmitoctove spektrum pasmoveho signalu v(t)')

%vykresleni fázového frekvenčního spektra
subplot(2,1,2);
stem(f_db,V_db_faze)
xlabel('f(Hz)')
ylabel('\Theta_m')
grid on
title('DSB-SC modulace - Fazove kmitoctove spektrum pasmoveho signalu v(t)')



