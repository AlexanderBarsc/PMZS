close all
clear all

syms t
%% Vykresleni signalu

amp = 0.2;
omega = 10*pi
fi = -pi/2
f = omega/(2*pi)
T = 1/f
tvek = 0:0.001:2*T;

w = amp*cos(omega*tvek + fi);

plot(tvek,w);
title('Vykresleni signalu w(t)');
xlabel('Cas t(s)');
ylabel('Amplituda');

%% Spektrum


m=-10:10;
wm=1/T*(int((amp*cos(omega*t + fi)*exp(-j*m*2*pi/T*t)),t,0,T))
wwm=double(wm); %lepsi zobrazeni hodnoty
A=abs(wwm);      % 
fi=angle(wwm); 
P=A.^2;
figure()
subplot(2,1,1)
stem(m*2*pi/T,A) %amplitudove spektrum
xlabel('\omega[rad/s]')
ylabel('|Wm|')
title('Amplitudove spektrum')
subplot(2,1,2)
stem(m*2*pi/T,fi) %fazove spektrum
xlabel('\omega[rad/s]')
ylabel('\phi')
title('Fazove spektrum')


%% AM modulace

f_c = 50;
R_t=abs(1+w); % amplitudova modulacni slozka 
if (w>-1)     % f√°zov√° modulaƒçn√≠ slo≈æka
    Theta_t = 0;
else 
    Theta_t = pi;
end

v_am_t=R_t.*cos(2*pi*f_c*tvek+Theta_t); 
figure()
hold on
xlabel('Cas t (s)');
ylabel('Amplituda');
title('AM modulace');
plot(tvek,v_am_t)
plot(tvek,w);
legend('Modulacni signal', 'Modulovany signal'); 

pocet_bodu = 3;
N_am = length(v_am_t);                                    
v_am_f = (fftshift(fft(v_am_t)))./N_am;                   
k_am = -N_am/2:N_am/2-1;                                 
f_am = k_am.*pocet_bodu*f_c./N_am;                        

V_am_amp = abs(v_am_f);                                   

for q=1:N_am
    if ((abs(real(v_am_f(q))) < 1e-3) && (abs(imag(v_am_f(q))) < 1e-3)) 
        v_am_f(q)=0;                                  
    end;
end;

V_am_faze = angle(v_am_f);

figure;
subplot(2,1,1);
stem(f_am,V_am_amp)
xlabel('f(Hz)')
ylabel('^F^R|W_m|')
grid on;
title('AM modulace - Amplitudova frekvencni spektrum pasmoveho signalu v(t)')

subplot(2,1,2);
stem(f_am,V_am_faze)
xlabel('f(Hz)')
ylabel('\Theta_m')
grid on
title('AM modulace - Fazove frekvencni spektrum pasmoveho signalu v(t)')


%% DSB-SC modulace - definov·nÌ a v˝poËet p·smovÈho sign·lu v(t)
% amplitudova modulacni slozka R_t = |m_t|
% fazova modulacni slozka Theta_t=0 pro m(t)>0, Theta_t=pi pro m(t)<0;
x_t = w; % souf·zov· modulacni slozka
y_t = 0;   % kvadraturnÌ modulacni slozka
% p·smov˝ (modulovan˝) sign·l v(t)
v_t = x_t.*cos(2*pi*f_c*tvek)-y_t.*sin(2*pi*f_c*tvek); 

%--------------------------------------------------------------------------
%V˝poËet frekvenËnÌho spektra p·smovÈho sign·lu pomocÌ funkce algoritmu FFT
N=length(v_t); % poËet hodnot p·smovÈho sign·lu v(t) 
v_f = (fftshift(fft(v_t)))./N; % komplexnÌ vektor frekvenËnÌho spektra v(t) 
k=-N/2:N/2-1; % pomocn˝ v˝poËet symetrickÈho pole osy x
f = k.*pocet_bodu*f_c./N; % v˝poËet x-ovÈ osy - frekvence f[Hz]

V_amp = abs(v_f);           %amplitudovÈ frekvenËnÌ spektrum v(t)
V_vykon = V_amp.^2;         %v˝konovÈ frekvenËnÌ spektrum v(t)

% Je-li nulov· hodnota amplitudy na danÈ frekvenci je takÈ f·ze = 0
for q=1:N
    if ((abs(real(v_f(q))) < 1e-3) && (abs(imag(v_f(q))) < 1e-3)) 
        v_f(q)=0;   
    end;
end;
V_faze = angle(v_f);        %f·zovÈ frekvenËnÌ spektrum v(t)

omega=2*pi*f;   % v˝poËet uhlovÈho kmitoËtu, omega=2*pi*f=2*pi/T


% GRAFY - ZOBRAZENÕ SIGN¡LŸ V »ASOV… OBLASTI
%--------------------------------------------------------------------------
% VykreslenÌ grafickÈho ËasovÈho pr˘bÏhu modulaËnÌho sign·lu
figure;
subplot(2,1,1);
plot(tvek,w);
title('ModulacnÌ sign·l w(t)');
ylabel('w(t)');
xlabel('t[s]');
 
%--------------------------------------------------------------------------
% VykreslenÌ grafickÈho ËasovÈho pr˘bÏhu p·smovÈho (modulovanÈho) sign·lu
subplot(2,1,2);
plot(tvek,v_t);
title('DSB-SC modulace - P·smov˝ sign·l v(t)');
ylabel('v(t)');
xlabel('t[s]');

%--------------------------------------------------------------------------
% GRAFY - ZOBRAZENÕ SIGN¡LŸ VE FREKVEN»NÕ OBLASTI
%--------------------------------------------------------------------------
%vykresleni amplitudovÈho frekvenËnÌho spektra 
figure;
subplot(3,1,1);
stem(f,V_amp)
xlabel('f[Hz]')
ylabel('^F^R|W_m|')
grid on;
title('DSB-SC modulace - AmplitudovÈ kmitoctove spektrum pasmoveho (modulovaneho)signalu v(t)')
%--------------------------------------------------------------------------
%vykresleni f·zovÈho frekvenËnÌho spektra
subplot(3,1,2);
stem(f,V_faze)
xlabel('f[Hz]')
ylabel('\Theta_m')
grid on
title('DSB-SC modulace - Fazove kmitoctove spektrum pasmoveho (modulovaneho)signalu v(t)')
%--------------------------------------------------------------------------
%vykresleni v˝konovÈho frekvenËnÌho spektra
subplot(3,1,3);
stem(f,V_vykon)
xlabel('f[Hz]')
ylabel('^F^R|P_m|') 
grid on
title('DSB-SC modulace - Vykonove kmitoctove spektrum pasmoveho (modulovaneho)signalu v(t)')
