clear all
close all

%% Vykresleni diskretniho signalu

N = 10; % Pocet vzorku
Fs = 25; % Vzorkovaci frekvence
Ts = 1/Fs; % Vzorkovaci perioda

Wmax = 0.2; % Amplituda
omega = 10*pi; % Uhlovy kmitocet
phase = -pi; % Fazovy posun

n = 0:1:(N - 1); % N potrebne pro jednu periodu
n2 = 0:1:(2*N - 1)

w = Wmax*cos(omega*n*Ts+phase) + 0.5
wdelsi = Wmax*cos(omega*n2*Ts+phase) + 0.5
stem(n,w);
xlabel('n');
ylabel('w[n]');
title('Harmonick� sign�l w[n]');

W_k1=fft(w);
W_k=fftshift(W_k1);

R_zaokrouhleno=round(real(W_k),6)
Im_zaokrouhleno=round(imag(W_k),6)

k = -(N/2):1:N/2-1;

figure;
subplot(2,1,1)
stem(k,(abs(W_k)/N))
xlabel('k[-]')
ylabel('|W[k]|')
title(' Amplitudov� frekven�n� spektrum |W[k]|')

subplot(2,1,2)
stem(k,(atan2(Im_zaokrouhleno,R_zaokrouhleno)))
xlabel('k[-]')
ylabel('phase W[k]')
title(' F�zov� frekven�n� spektrum W[k]')

%% Obdelnikova okenni funkce

[Wo, no] = VygenerujObdelnikovouOkenniFunkci(length(n), 0, N - 1)
Wres = Wo.*w

figure()
stem(no, Wres)
title('Signal w[n] vynasobeny s obdelnikovou okenni funkci')
xlabel('n');
ylabel('wObd�ln�k[n]');


W_kobdelnik = fftshift(fft(Wres));

R_zaokrouhleno=round(real(W_kobdelnik),6)
Im_zaokrouhleno=round(imag(W_kobdelnik),6)

figure;
subplot(2,1,1)
stem(k,(abs(W_kobdelnik)/length(no)))
xlabel('k[-]')
ylabel('|W[k]|')
title(' Amplitudov� frekven�n� spektrum sign�lu vyn�soben�ho s obdeln�kovou okenn� funkc� |W[k]|')

subplot(2,1,2)
stem(k,(atan2(Im_zaokrouhleno,R_zaokrouhleno)))
xlabel('k[-]')
ylabel('phase W[k]')
title(' F�zov� frekven�n� spektrum sign�lu vyn�soben�ho s obdeln�kovou okenn� funkc� W[k]')


%% Hanningova okenn� funkce

nHanning = 0:1:N-1;
WHanning = 1 - cos((2*pi*nHanning)/(N - 1))

Wres = WHanning .* w;

figure()
stem(n, Wres)

title('Sign�l w[n] vyn�soben� s Hanningovou okenn� funkc�');
xlabel('n');
ylabel('wHanning[n]');

W_kHanning = fftshift(fft(Wres));

R_zaokrouhleno=round(real(W_kHanning),6)
Im_zaokrouhleno=round(imag(W_kHanning),6)

figure;
subplot(2,1,1)
stem(k,(abs(W_kHanning)/length(nHanning)))
xlabel('k[-]')
ylabel('|W[k]|')
title(' Amplitudov� frekven�n� spektrum sign�lu vyn�soben�ho s Hanningovou okenn� funkc� |W[k]|')

subplot(2,1,2)
stem(k,(atan2(Im_zaokrouhleno,R_zaokrouhleno)))
xlabel('k[-]')
ylabel('phase W[k]')
title('F�zov� frekven�n� spektrum sign�lu vyn�soben�ho s Hanningovou okenn� funkc� W[k]')


%% Flat top okenni funkce

nFlatTop = 0:1:N-1;
WFlatTop = 1 - 1.98*cos((2*pi*nFlatTop)/(N - 1)) + 1.29*cos((4*pi*nFlatTop)/(N - 1)) - 0.388*cos((6*pi*nFlatTop)/(N - 1)) + 0.0322*cos((8*pi*nFlatTop)/(N - 1))

Wres = WFlatTop .* w;

figure()
stem(n, Wres)

title('Sign�l w[n] vyn�soben� s FlatTop okenn� funkc�');
xlabel('n');
ylabel('wFlatTop[n]');

W_kFlatTop = fftshift(fft(Wres));

R_zaokrouhleno=round(real(W_kFlatTop),6)
Im_zaokrouhleno=round(imag(W_kFlatTop),6)

figure;
subplot(2,1,1)
stem(k,(abs(W_kFlatTop)/length(nFlatTop)))
xlabel('k[-]')
ylabel('|W[k]|')
title(' Amplitudov� frekven�n� spektrum sign�lu vyn�soben�ho s FlatTop okenn� funkc� |W[k]|')

subplot(2,1,2)
stem(k,(atan2(Im_zaokrouhleno,R_zaokrouhleno)))
xlabel('k[-]')
ylabel('phase W[k]')
title('F�zov� frekven�n� spektrum sign�lu vyn�soben�ho s FlatTop okenn� funkc� W[k]')

%% Barlettova okenn� funkce

nBarlett = 0:1:N-1;
WBarlett = zeros(1, length(nBarlett))

for i = 0:length(nBarlett) - 1
   
    if(i <= (N - 1)/2)
        WBarlett(1,i + 1) = 2*i/(N - 1);
    else
       WBarlett(1,i + 1) = 2 - (2*i/(N - 1));
    end
    
end

Wres = WBarlett .* w;

figure()
stem(n, Wres)

title('Sign�l w[n] vyn�soben� s Barlettovou okenn� funkc�');
xlabel('n');
ylabel('wBarlett[n]');

W_kBarlett = fftshift(fft(Wres));

R_zaokrouhleno=round(real(W_kBarlett),6)
Im_zaokrouhleno=round(imag(W_kBarlett),6)

figure;
subplot(2,1,1)
stem(k,(abs(W_kBarlett)/length(nBarlett)))
xlabel('k[-]')
ylabel('|W[k]|')
title(' Amplitudov� frekven�n� spektrum sign�lu vyn�soben�ho s Barlettovou okenn� funkc� |W[k]|')

subplot(2,1,2)
stem(k,(atan2(Im_zaokrouhleno,R_zaokrouhleno)))
xlabel('k[-]')
ylabel('phase W[k]')
title('F�zov� frekven�n� spektrum sign�lu vyn�soben�ho s Barlettovou okenn� funkc� W[k]')

%% Hannova okenn� funkce

nHannova = 0:1:N -1;
WHannova = 0.5 - 0.5*cos((2*pi*nHannova)/(N - 1))

Wres = WHannova .* w;

figure()
stem(n, Wres)

title('Sign�l w[n] vyn�soben� s Hannovou okenn� funkc�');
xlabel('n');
ylabel('wHannova[n]');

W_kHannova = fftshift(fft(Wres));

R_zaokrouhleno=round(real(W_kHannova),6)
Im_zaokrouhleno=round(imag(W_kHannova),6)

figure;
subplot(2,1,1)
stem(k,(abs(W_kHannova)/length(nHannova)))
xlabel('k[-]')
ylabel('|W[k]|')
title(' Amplitudov� frekven�n� spektrum sign�lu vyn�soben�ho s Hannovou okenn� funkc� |W[k]|')

subplot(2,1,2)
stem(k,(atan2(Im_zaokrouhleno,R_zaokrouhleno)))
xlabel('k[-]')
ylabel('phase W[k]')
title('F�zov� frekven�n� spektrum sign�lu vyn�soben�ho s Hannovou okenn� funkc� W[k]')






