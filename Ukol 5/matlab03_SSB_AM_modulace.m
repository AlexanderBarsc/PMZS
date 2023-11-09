%--------------------------------------------------------------------------
%           Výuková podpora pøemìtu Modulované signály
%--------------------------------------------------------------------------
%                       Program SSB-AM MODULACE (USSB-AM, LSSB-AM)
%--------------------------------------------------------------------------
% 
% Výpoèet a grafické vykreslení modulovaného vysokofrekvenèního signálu 
% pomocí analogové amplitudové SSB-AM
% (SingleSide Band - Amplitude Modulation)
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
pocet_bodu=3; %poèet bodù v èase bìhem jedné periody fc
deltat=1/(f_c*pocet_bodu); % zanedbatelný èasový okamžik
t=0:deltat:cas-deltat; % èasový prùbìh

%--------------------------------------------------------------------------
% Definování a výpoèet modulaèního signálu m(t)
% m_t je složeny ze dvou harmonických složek

Wmax1 = 3;    % amplituda 1.harmonické složky modulaèního signálu m(t)
f_m1 = 100;   % frekvence 1.harmonické složky modulaèního signálu m(t)
faze1 = pi/8; % fáze 1.harmonické složky modulaèního signálu m(t)
Wmax2 = 5;    % amplituda 2.harmonické složky modulaèního signálu m(t)
f_m2 = 40;    % frekvence 2.harmonické složky modulaèního signálu m(t)
faze2 = pi;   % fáze 2.harmonické složky modulaèního signálu m(t)
% informaèní (modulaèní) signál m(t)
m_t=Wmax1*cos(2*pi*f_m1*t+faze1)+Wmax2*cos(2*pi*f_m2*t+faze2); 
%--------------------------------------------------------------------------
% výpoèet Hilbertovy transformace signalu m(t)
m_H_t=Wmax1*sin(2*pi*f_m1*t+faze1)+Wmax2*sin(2*pi*f_m2*t+faze2); 

%--------------------------------------------------------------------------
% SSB-AM modulace - definování a výpoèet pásmového signálu v(t)
% komplexní obálka g(t) = m(t)+-jmH(t)
% v(t)= Re(g(t)).exp(j.omega.fc.t)

% 1.USSB-AM modulace - Upper Single Side Band - Amplitude Modulation
%v_t = real((m_t+j*m_H_t).*exp(j*2*pi*f_c*t)); % pásmový (modulovaný) signál v(t)

% 2.LSSB-AM modulace - Lower Single Side Band - Amplitude Modulation
v_t = real((m_t-j*m_H_t).*exp(j*2*pi*f_c*t)); % pásmový (modulovaný) signál v(t)

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
subplot(3,1,1);
plot(t,m_t);
title('Modulacní signál m(t)');
ylabel('m(t)');
xlabel('t[s]');
 
%--------------------------------------------------------------------------
% Vykreslení grafického èasového prùbìhu transformovaného modulaèního signálu
% po Hilbertovì transformaci
subplot(3,1,2);
plot(t,m_H_t);
title('Transformovaný modulacní signál po Hilbertove transformaci m_H(t)');
ylabel('m_H(t)');
xlabel('t[s]');

%--------------------------------------------------------------------------
% Vykreslení grafického èasového prùbìhu pásmového (modulovaného) signálu
subplot(3,1,3);
plot(t,v_t);
title('SSB-AM modulace - Pásmový signál v(t)');
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
title('SSB-AM modulace - Amplitudové frekvencni spektrum pasmoveho (modulovaneho)signalu v(t)')
%--------------------------------------------------------------------------
%vykresleni fázového frekvenèního spektra
subplot(3,1,2);
stem(f,V_faze)
xlabel('f[Hz]')
ylabel('\Theta_m')
grid on
title('SSB-AM modulace - Fazove frekvencni spektrum pasmoveho (modulovaneho)signalu v(t)')
%--------------------------------------------------------------------------
%vykresleni výkonového frekvenèního spektra
subplot(3,1,3);
stem(f,V_vykon)
xlabel('f[Hz]')
ylabel('^F^R|P_m|') 
grid on
title('SSB-AM modulace - Vykonove frekvencni spektrum pasmoveho (modulovaneho)signalu v(t)')

