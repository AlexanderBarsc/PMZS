clear all
close all


%% Cviceni 3 vypracovani

% Vykreslene signaly m1 a m2

tsampling = 0.001

t = 0:tsampling:1
f1 = 5;
omega1 = 2*pi*f1;
U1 = 2;

f2 = 8;
omega2 = 2*pi*f2;
U2 = 3;

m1 = U1*cos(omega1*t);

m2 = U2*cos(omega2*t);

figure();
plot(t,m1);
hold on
plot(t,m2);
legend('m1', 'm2')
ylabel('Amplituda');
xlabel('Cas (t)');
title('Vykreslene signaly m1 a m2');
xlim([0 0.250])

%% Soucet signalu m1 a m2 -> m3

m3 = m1 + m2;

tmaxharmonic = 1/f2;
tsamplingNew = tmaxharmonic/5; % aspon pet vzorku na periodu
tfit = tsamplingNew / tsampling; %

tMaxHarmonic2periods = (2*tmaxharmonic)/tsampling; % 2 periody nejkratsi periody harmonicky slozky


m3new = m3(1:tfit:251); % Vytahneme kazdy petinovy vzorek z periody

t3new = (0:tfit:tfit*10)/1000; 


figure()
plot(t, m3)
ylabel('Amplituda');
xlabel('Cas (t)');
title('Soucet signalu m1 a m2 = signal m3');

figure()
stem(t3new, m3new)
hold on
plot(t, m3)
ylabel('Amplituda');
xlabel('Cas (t)');
title('Puvodni signal a navzorkovany signal m3');
legend('Navzorkovany m3', 'Puvodni m3')
xlim([0 0.250])




%% PCM kvantovani

figure()
hold on
krok = 10/15;
pcmhladiny = -5:krok:5;

for j = 1: length(pcmhladiny)
   
    yline(pcmhladiny(j));
    
end



for i = 1: length(m3new)
    
    for j = 1: length(pcmhladiny)
        
        if(j == length(pcmhladiny))
            
          pcmkvant(i) = pcmhladiny(j)
            break;
        end 
        
        if(m3new(i) >= pcmhladiny(j) && m3new(i) < pcmhladiny(j + 1))
            
            pcmkvant(i) = pcmhladiny(j)
            break;
        end
    
    end
    
end

ylabel('Amplituda');
xlabel('Cas (t)');
title('Kvantovany signal zobrazeny s originalnim signalem a hladinami');
plot(t,m3);
stem(t3new, pcmkvant)
xlim([0 0.250])
hold on

%% PCM vystup 

t3pcm = zeros(1,2*length(m3new))
m3pcm = zeros(1,2*length(m3new))

i = 1
j = 1;
while i <= length(pcmkvant)
    
    if(i == length(pcmkvant))
    t3pcm(j) = t3new(i);
    m3pcm(j) = pcmkvant(i);
    t3pcm(j + 1) = t3new(i);
    m3pcm(j + 1) = pcmkvant(i);
    break; 
        
    end
    
    t3pcm(j) = t3new(i);
    m3pcm(j) = pcmkvant(i);
    t3pcm(j + 1) = t3new(i + 1);
    m3pcm(j + 1) = pcmkvant(i);
    
    i = i + 1;
    j = j + 2;

end
    

figure()
hold on
plot(t,m3);
stem(t3new, pcmkvant)
plot(t3pcm, m3pcm);
xlim([0 0.250])
ylabel('Amplituda');
xlabel('Cas (t)');
title('PCM vystup zobrazeny s puvodnim signalem a kvantovanymi vzorky');
legend('Puvodni signal', 'Kvantovany signal', 'PCM vystup');
%% PAM modulace

m3pam(1) = m3new(1)
t3pam(1) = 0

for i = 1:length(m3new)

    for j = 1:25

        if( j < 10)
            m3pam(end + 1) = m3new(i);
        elseif(j == 10)
            m3pam(end + 1) = m3new(i);
            t3pam(end + 1) = t3pam(end) + 0.001;
            m3pam(end + 1) = 0;
            t3pam(end + 1) = t3pam(end);
            continue;
        elseif(j == 25)

            if(i == length(m3new))
               break; 
            end
            m3pam(end + 1) = 0;
            t3pam(end + 1) = t3pam(end) + 0.001;
            m3pam(end + 1) = m3new(i + 1);
            t3pam(end + 1) = t3pam(end);
            break;
        else
            m3pam(end + 1) = 0;

        end;

            t3pam(end + 1) = t3pam(end) + 0.001;

    end

end

figure()
plot(t3pam, m3pam)
hold on
stem(t3new, m3new)
plot(t,m3);
xlim([0 0.250])

legend('PAM modulace', 'm3 vzorkovany', 'Puvodni signal m3')
ylabel('Amplituda');
xlabel('Cas (t)');
title('PAM modulace aplikovana na vzorkovany signal m3');



