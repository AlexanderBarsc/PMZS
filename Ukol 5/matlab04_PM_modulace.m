%--------------------------------------------------------------------------
%           Výuková podpora pøemìtu Modulované signály
%--------------------------------------------------------------------------
%                       Program PM MODULACE 
%--------------------------------------------------------------------------
% 
% Výpoèet a grafické vykreslení modulovaného vysokofrekvenèního signálu 
% pomocí analogové úhlové PM modulace (Phase Modulation)
%
% Volitelný modulaèní (informaèní) signál a vysokofrekvenèní nosná
% 
% VŠB - Technická univerzita Ostrava
% Fakulta elektrotechniky a informatiky
% Katedra mìøicí a øídicí techniky
% 17.listopadu 15
% Ostrava - Poruba
% 708 33
%
% Vypracoval Zdenìk Macháèek 2010
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Vyèištìní všech promìnných pamìti a uzavøení všech grafických oken
clear all;
close all;

%--------------------------------------------------------------------------
% Specifikace èasového prùbìhu-èas t[s] a nosné frekvence fc[Hz]RF signálu
f_c=200; % nosná frekvence vysokofrekvenèního RF signálu
cas = 0.2; % èas prùbìhu simulace 
pocet_bodu=10; %poèet bodù v èase bìhem jedné periody fc
deltat=1/(f_c*pocet_bodu); % zanedbatelný èasový okamžik
t=0:deltat:cas-deltat; % èasový prùbìh

%--------------------------------------------------------------------------
% Definování a výpoèet modulaèního signálu m(t)
% m_t je volitelný harmonický nebo obdelníkový signál

%{
%1.Hamonický signál m(t)
Wmax = 1.5;  % amplituda modulaèního signálu m(t)
f_m = 10;    % frekvence modulaèního signálu m(t)
faze = pi/2; % fáze modulaèního signálu m(t)
m_t=Wmax*cos(2*pi*f_m*t+faze); % informaèní (modulaèní) signál m(t)
%}

%2.Obdelníkový signál m(t)
Wmax = 1.5;  % velikost modulaèního obdelníkového signálu m(t)
t_p1 =length(t)-300; % poèet prvkù v èase pro hodnotu m1(t) = 0
m1_t(1:t_p1)=0;      % první èást signálu m1(t)
t_p2 =length(t)-200; % poèet prvkù v èase pro hodnotu m2(t) = Wmax
m2_t(1:t_p2)=Wmax;   % druhá èást signálu m2(t)
t_p3 =((length(t)-300)); % poèet prvkù v èase pro hodnotu m3(t) = 0
m3_t(1:t_p3)=0;      % tøetí èást signálu m3(t)
m_t=[m1_t m2_t m3_t]; % výsledný složený modulaèní signál m(t)


%--------------------------------------------------------------------------
% PM modulace - definování a výpoèet pásmového signálu v(t)
% amplitudova modulacni slozka R_t = 1
% fazova modulacni slozka Theta_t=Dp.m(t)

A_c = 2;            % zesílení amplitudy modulovaného signálu v(t)
R_t=1;              % amplitudova modulacni slozka 
Dp = 1;             % Index fázové modulace
Theta_t = Dp*m_t;   % Fázová modulaèní složka
omega_c = 2*pi*f_c; % uhlový kmitoèet nosného signálu
v_t = A_c*R_t*cos(omega_c*t+Theta_t); % pásmový (modulovaný) signál v(t)

%--------------------------------------------------------------------------
%Výpoèet frekvenèního spektra pásmového signálu pomocí funkce algoritmu FFT
N=length(v_t); % poèet hodnot pásmového signálu v(t) 
v_f = (fftshift(fft(v_t)))./N; % komplexní vektor frekvenèního spektra v(t) 
k=-N/2:N/2-1; % pomocný výpoèet symetrického pole osy x
f = k.*pocet_bodu*f_c./N; % výpoèet x-ové osy - frekvence f[Hz]

V_amp = abs(v_f);           %amplitudové frekvenèní spektrum v(t)
V_vykon = V_amp.^2;         %výkonové frekvenèní spektrum v(t)

% Je-li nulová hodnota amplitudy na dané frekvenci je také fáze = 0
for q=1:N
    if ((abs(real(v_f(q))) < 3e-2) && (abs(imag(v_f(q))) < 3e-2)) 
        v_f(q)=0;   
    end;
end;
V_faze = angle(v_f);        %fázové frekvenèní spektrum v(t)

%--------------------------------------------------------------------------
% GRAFY - ZOBRAZENÍ SIGNÁLÙ V ÈASOVÉ OBLASTI
%--------------------------------------------------------------------------
% Vykreslení grafického èasového prùbìhu modulaèního signálu
figure;
subplot(2,1,1);
plot(t,m_t);
title('Modulacní signál m(t)');
ylabel('m(t)');
xlabel('t[s]');
 
%--------------------------------------------------------------------------
% Vykreslení grafického èasového prùbìhu pásmového (modulovaného) signálu
subplot(2,1,2);
plot(t,v_t);
title('PM modulace - Pásmový signál v(t)');
ylabel('v(t)');
xlabel('t[s]');

%--------------------------------------------------------------------------
% GRAFY - ZOBRAZENÍ SIGNÁLÙ VE FREKVENÈNÍ OBLASTI
%--------------------------------------------------------------------------
%vykresleni amplitudového frekvenèního spektra 
figure;
subplot(3,1,1);
stem(f,V_amp)
xlabel('f[Hz]')
ylabel('^F^R|W_m|')
grid on;
title('PM modulace - Amplitudové frekvencni spektrum pasmoveho (modulovaneho)signalu v(t)')
%--------------------------------------------------------------------------
%vykresleni fázového frekvenèního spektra
subplot(3,1,2);
stem(f,V_faze)
xlabel('f[Hz]')
ylabel('\Theta_m')
grid on
title('PM modulace - Fazove frekvencni spektrum pasmoveho (modulovaneho)signalu v(t)')
%--------------------------------------------------------------------------
%vykresleni výkonového frekvenèního spektra
subplot(3,1,3);
stem(f,V_vykon)
xlabel('f[Hz]')
ylabel('^F^R|P_m|') 
grid on
title('PM modulace - Vykonove frekvencni spektrum pasmoveho (modulovaneho)signalu v(t)')
