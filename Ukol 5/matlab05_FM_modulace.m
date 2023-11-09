%--------------------------------------------------------------------------
%           Výuková podpora pøemìtu Modulované signály
%--------------------------------------------------------------------------
%                       Program FM MODULACE 
%--------------------------------------------------------------------------
% 
% Výpoèet a grafické vykreslení modulovaného vysokofrekvenèního signálu 
% pomocí analogové úhlové FM modulace (Frequency Modulation)
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
f_c=100; % nosná frekvence vysokofrekvenèního RF signálu
cas = 0.4; % èas prùbìhu simulace 
pocet_bodu=10; %poèet bodù v èase bìhem jedné periody fc
deltat=1/(f_c*pocet_bodu); % zanedbatelný èasový okamžik
t=0:deltat:cas-deltat; % èasový prùbìh

%--------------------------------------------------------------------------
% Definování a výpoèet modulaèního signálu m(t)
% m_t je volitelný harmonický nebo obdelníkový signál

%{
%1.Hamonický signál m(t)
Wmax = 4;  % amplituda modulaèního signálu m(t)
f_m = 2;    % frekvence modulaèního signálu m(t)
m_t=Wmax*cos(2*pi*f_m*t); % informaèní (modulaèní) signál m(t)
%}

%%{
%2.Obdelníkový signál m(t)
Wmax = 1.5;  % velikost modulaèního obdelníkového signálu m(t)
t_p1 =length(t)-300; % poèet prvkù v èase pro hodnotu m1(t) = 0
m1_t(1:t_p1)=0;      % první èást signálu m1(t)
t_p2 =length(t)-200; % poèet prvkù v èase pro hodnotu m2(t) = Wmax
m2_t(1:t_p2)=Wmax;   % druhá èást signálu m2(t)
t_p3 =((length(t)-300)); % poèet prvkù v èase pro hodnotu m3(t) = 0
m3_t(1:t_p3)=0;      % tøetí èást signálu m3(t)
m_t=[m1_t m2_t m3_t]; % výsledný složený modulaèní signál m(t)
%}

%--------------------------------------------------------------------------
% FM modulace - definování a výpoèet pásmového signálu v(t)
% amplitudova modulacni slozka R_t = 1
% fazova modulacni slozka Theta_t=Df.integrál{m(t)}

A_c = 2;     % zesílení amplitudy modulovaného signálu v(t)
R_t = 1;     % amplitudova modulacni slozka 
%Df = 50;     % 1.Index fázové modulace pro harmonický signál m(t)
Df = 200;   % 2.Index fázové modulace pro obdelníkový signál m(t)
omega_c = 2*pi*f_c; % uhlový kmitoèet nosného signálu

%1.Výpoèet integrálu pro harmonický signál m(t)
%Theta_t = Df*(Wmax/(2*pi*f_m)*sin((2*pi*f_m*t))); % Fázová modulaèní složka

%%{
%2.Výpoèet integrálu pro obdelníkový signál m(t)
Theta_t1(1:t_p1) = 0;     % první èást fázové modulaèní složky
t2 = 0:deltat:(t_p2*deltat)-deltat;
Theta_t2(1:t_p2) = Df*t2; % druhá èást fázové modulaèní složky
Theta_t3(1:t_p3) = Df*t_p2;     % tøetí èást fázové modulaèní složky
Theta_t=[Theta_t1 Theta_t2 Theta_t3]; % výsledná složená fázová modulaèní složka Theta(t)
%}

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
    if ((abs(real(v_f(q))) < 1e-2) && (abs(imag(v_f(q))) < 1e-2)) 
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
title('FM modulace - Pásmový signál v(t)');
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
title('FM modulace - Amplitudové frekvencni spektrum pasmoveho (modulovaneho)signalu v(t)')
%--------------------------------------------------------------------------
%vykresleni fázového frekvenèního spektra
subplot(3,1,2);
stem(f,V_faze)
xlabel('f[Hz]')
ylabel('\Theta_m')
grid on
title('FM modulace - Fazove frekvencni spektrum pasmoveho (modulovaneho)signalu v(t)')
%--------------------------------------------------------------------------
%vykresleni výkonového frekvenèního spektra
subplot(3,1,3);
stem(f,V_vykon)
xlabel('f[Hz]')
ylabel('^F^R|P_m|') 
grid on
title('FM modulace - Vykonove frekvencni spektrum pasmoveho (modulovaneho)signalu v(t)')
