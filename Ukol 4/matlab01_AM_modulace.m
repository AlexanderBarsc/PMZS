%--------------------------------------------------------------------------
%           Výuková podpora pøemìtu Modulované signály
%--------------------------------------------------------------------------
%                       Program AM MODULACE 
%--------------------------------------------------------------------------
% 
% Výpoèet a grafické vykreslení modulovaného vysokofrekvenèního signálu 
% pomocí analogové amplitudové AM modulace (Amplitude Modulation)
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
f_c=1000; % nosná frekvence vysokofrekvenèního RF signálu
cas = 0.1; % èas prùbìhu simulace 
pocet_bodu=30; %poèet bodù v èase bìhem jedné periody fc
deltat=1/(f_c*pocet_bodu); % zanedbatelný èasový okamžik
t=0:deltat:cas-deltat; % èasový prùbìh

%--------------------------------------------------------------------------
% Definování a výpoèet modulaèního signálu m(t)
% m_t resp. Wmax musí být menší než 1, jinak je signál v(t) pøemodulovaný


Wmax1 = 0.5;  % amplituda 1.harmonické složky modulaèního signálu m(t)
f_m1 = 10;    % frekvence 1.harmonické složky modulaèního signálu m(t)
faze1 = pi/2; % fáze 1.harmonické složky modulaèního signálu m(t)
Wmax2 = 0.3;  % amplituda 2.harmonické složky modulaèního signálu m(t)
f_m2 = 20;    % frekvence 2.harmonické složky modulaèního signálu m(t)
faze2 = pi/4; % fáze 2.harmonické složky modulaèního signálu m(t)
m_t=Wmax1*cos(2*pi*f_m1*t+faze1)+Wmax2*cos(2*pi*f_m2*t+faze2); % informaèní (modulaèní) signál m(t)

%{
Wmax = 0.75;  % amplituda modulaèního signálu m(t)
f_m = 40;    % frekvence modulaèního signálu m(t)
faze = pi/2; % fáze modulaèního signálu m(t)
m_t=Wmax*cos(2*pi*f_m*t+faze); % informaèní (modulaèní) signál m(t)

%}

%--------------------------------------------------------------------------
% AM modulace - definování a výpoèet pásmového signálu v(t)
% amplitudova modulacni slozka R_t = |1+m_t|
% fazova modulacni slozka Theta_t=0 pro m(t)>-1, Theta_t=pi pro m(t)<-1;
% m(t)nesmi byt premodulovany >1

R_t=abs(1+m_t); % amplitudova modulacni slozka 
if (m_t>-1)     % fázová modulaèní složka
    Theta_t = 0;
else 
    Theta_t = pi;
end
v_t=R_t.*cos(2*pi*f_c*t+Theta_t); % pásmový (modulovaný) signál v(t)

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
title('AM modulace - Pásmový signál v(t)');
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
title('AM modulace - Amplitudové frekvencni spektrum pasmoveho (modulovaneho)signalu v(t)')
%--------------------------------------------------------------------------
%vykresleni fázového frekvenèního spektra
subplot(3,1,2);
stem(f,V_faze)
xlabel('f[Hz]')
ylabel('\Theta_m')
grid on
title('AM modulace - Fazove frekvencni spektrum pasmoveho (modulovaneho)signalu v(t)')
%--------------------------------------------------------------------------
%vykresleni výkonového frekvenèního spektra
subplot(3,1,3);
stem(f,V_vykon)
xlabel('f[Hz]')
ylabel('^F^R|P_m|') 
grid on
title('AM modulace - Vykonove frekvencni spektrum pasmoveho (modulovaneho)signalu v(t)')












