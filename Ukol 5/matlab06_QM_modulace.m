%--------------------------------------------------------------------------
%           Výuková podpora pøemìtu Modulované signály
%--------------------------------------------------------------------------
%                       Program QM MODULACE 
%--------------------------------------------------------------------------
% 
% Výpoèet a grafické vykreslení modulovaného vysokofrekvenèního signálu 
% pomocí analogové kvadraturní QM modulace (Quadrature Modulation)
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
f_c=5000; % nosná frekvence vysokofrekvenèního RF signálu
cas = 0.05; % èas prùbìhu simulace 
pocet_bodu=10; %poèet bodù v èase bìhem jedné periody fc
deltat=1/(f_c*pocet_bodu); % zanedbatelný èasový okamžik
t=0:deltat:cas-deltat; % èasový prùbìh

%--------------------------------------------------------------------------
% Definování a výpoèet modulaèního signálu m(t)
% m_t je volitelný harmonický nebo obdelníkový signál

%1.Hamonický signál m(t)
Wmax1 = 4;  % amplituda modulaèního signálu m(t)
f_m1 = 100;    % frekvence modulaèního signálu m(t)
m_t1 = Wmax1*cos(2*pi*f_m1*t); % informaèní (modulaèní) signál m(t)

%2.Hamonický signál m(t)
Wmax2 = 6;  % amplituda modulaèního signálu m(t)
f_m2 = 400;    % frekvence modulaèního signálu m(t)
m_t2 = Wmax2*cos(2*pi*f_m2*t); % informaèní (modulaèní) signál m(t)


%--------------------------------------------------------------------------
% QM modulace - definování a výpoèet pásmového signálu v(t)
omega_c = 2*pi*f_c; % uhlový kmitoèet nosného signálu

%Výpoèet komplexní obálky signálu v(t)=m1(t)+j.m2(t)
g_t = m_t1 + 1j*m_t2;          % komplexni obalka

%Vypocet pasmoveho (modulovaneho) signalu z algoritmu komplexní obálky
v_t = real(g_t.*exp(1j*omega_c*t));

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
    if ((abs(real(v_f(q))) < 1e-3) && (abs(imag(v_f(q))) < 1e-3)) 
        v_f(q)=0;   
    end;
end;
V_faze = angle(v_f);        %fázové frekvenèní spektrum v(t)

%--------------------------------------------------------------------------
% GRAFY - ZOBRAZENÍ SIGNÁLÙ V ÈASOVÉ OBLASTI
%--------------------------------------------------------------------------
% Vykreslení grafického èasového prùbìhu modulaèního signálu m1(t)
figure;
subplot(3,1,1);
plot(t,m_t1);
title('Modulacní signál m_1(t)');
ylabel('m_1(t)');
xlabel('t[s]');
 
%--------------------------------------------------------------------------
% Vykreslení grafického èasového prùbìhu modulaèního signálu m2(t)
subplot(3,1,2);
plot(t,m_t2);
title('Modulacní signál m_2(t)');
ylabel('m_2(t)');
xlabel('t[s]');

%--------------------------------------------------------------------------
% Vykreslení grafického èasového prùbìhu pásmového (modulovaného) signálu
subplot(3,1,3);
plot(t,v_t);
title('QM modulace - Pásmový signál v(t)');
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
title('QM modulace - Amplitudové frekvencni spektrum pasmoveho (modulovaneho)signalu v(t)')
%--------------------------------------------------------------------------
%vykresleni fázového frekvenèního spektra
subplot(3,1,2);
stem(f,V_faze)
xlabel('f[Hz]')
ylabel('\Theta_m')
grid on
title('QM modulace - Fazove frekvencni spektrum pasmoveho (modulovaneho)signalu v(t)')
%--------------------------------------------------------------------------
%vykresleni výkonového frekvenèního spektra
subplot(3,1,3);
stem(f,V_vykon)
xlabel('f[Hz]')
ylabel('^F^R|P_m|') 
grid on
title('QM modulace - Vykonove frekvencni spektrum pasmoveho (modulovaneho)signalu v(t)')



