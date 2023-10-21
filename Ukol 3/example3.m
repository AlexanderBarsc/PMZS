
clear all;
close all;
clc;

%% --------------------------------------------------------------------------

U1 = 2;
Omega1 = 6*pi;
phi1 = pi;
T1 = 2*pi/Omega1;
t1 = 0:0.001:2*T1;

w1 = U1*cos(Omega1*t1+phi1);

figure
plot(t1, w1)
xlabel('t')
ylabel('w1(t)')
title('Harmonic signal w1(t)')
xlim([0, 2*T1])


U2 = 4;
Omega2 = 2*pi;
phi2 = -pi/2;
T2 = 2*pi/Omega2;
t2 = 0:0.001:2*T2;

w2 = U2*cos(Omega2*t2+phi2);

figure
plot(t2, w2)
xlabel('t')
ylabel('w2(t)')
title('Harmonic signal w2(t)')
xlim([0, 2*T2])

%% seceteni signalu



T = 1;
t = 0:0.001:2*T;



w3 = U1*cos(Omega1*t+phi1) + U2*cos(Omega2*t+phi2);

figure
plot(t, w3)
xlabel('t')
ylabel('w3(t)')
title('Signal w3(t)')
xlim([0, 2*T])


%% vzorkovani


pocet_vzorku_na_periodu = 10;
pocet_vykreslenych_period = 6;
Ts = 2*pi/(Omega1*pocet_vzorku_na_periodu);
N = pocet_vzorku_na_periodu * pocet_vykreslenych_period;
n = 0:1:N-1;

w3_n = U1*cos(Omega1*n*Ts+phi1) + U2*cos(Omega2*n*Ts+phi2);

figure
stem(n,w3_n); 
xlabel('n')
ylabel('w3_n[n]')
title('Signal w3_n[n]')


%% kvantovani

% pro cele cislo staci floor

w_kvant = floor(w3_n);

figure
stem(n,w_kvant); 
xlabel('n')
ylabel('w-kvant[n]')
title('Signal w-kvant[n]')


%% kodovani

pocet_bitu = 4;
offset = min(w_kvant);

w_kvant_offset = w_kvant - offset;

w_bin = zeros(1,(length(w_kvant_offset)*pocet_bitu));

for i = 1:(length(w_kvant_offset))
    current_number = w_kvant_offset(i);
    for j = 1:1:pocet_bitu
        binary = 0;
        index = pocet_bitu - j;
        if current_number >= 2^index
            binary = 1;
            current_number = current_number - 2^index;
        end
        w_bin((i-1)*pocet_bitu+j) = binary;
    end
end
n_bin = 0:1/pocet_bitu:N-1/pocet_bitu;

T = 1;
tick =0.001;
t = 0:tick:(2*T-tick)

w_kodovany = zeros(1, length(t));
for i = 1:1:length(t)
    index = t(i)*pocet_bitu/(Ts);
    w_kodovany(i) = w_bin(int32(floor(index))+1);
end


figure
plot(t, w_kodovany);
ylim([-1,2])
grid on;
xlabel('t')
ylabel('w-kodovany(t)')
title('Signal w-PCM(t)')


%% modulace pwm



offset = min(w_kvant);
maximum = max(w_kvant);
max_steps = maximum - offset;
w_kvant_offset = w_kvant - offset;

T = 1;
tick =0.001;
t = 0:tick:(2*T-tick);

w_kodovany = zeros(1, length(t));
for i = 1:1:length(t)
    index = t(i)/Ts;
    int_index = int32(floor(index));
    remainer = index - double(int_index);
    current_value = w_kvant_offset(int_index+1);
    percentage_requirement = current_value / max_steps;
    w_kodovany(i) = 0;
    if remainer < percentage_requirement
        w_kodovany(i) = 1;
    end
   

end


figure
plot(t, w_kodovany);
ylim([-1,2])
grid on;
xlabel('t')
ylabel('w-modulovany(t)')
title('Signal w-PWM(t)')



%% pnm

% nefunguje

offset = min(w_kvant);
maximum = max(w_kvant);
max_steps = maximum - offset;
w_kvant_offset = w_kvant - offset;

T = 1;
tick =0.001;
t = 0:tick:(2*T-tick);

impulse_width = 0.05;

w_kodovany = zeros(1, length(t));
for i = 1:1:length(t)
    index = t(i)/Ts;
    int_index = int32(floor(index));
    remainer = index - double(int_index);
    current_value = w_kvant_offset(int_index+1);
    percentage_requirement = current_value / max_steps;
    w_kodovany(i) = 0;
    if remainer < percentage_requirement
        pulse_polarity = t(i)/impulse_width;
        pulse_polarity_floor = floor(pulse_polarity);
        pulse_polarity_x = pulse_polarity - pulse_polarity_floor;
        if pulse_polarity_x < 0.5
            w_kodovany(i) = 1;
        end
        w_kodovany(i) = 1;
    end
   

end


figure
plot(t, w_kodovany);
ylim([-1,2])
grid on;
xlabel('t')
ylabel('w-modulovany(t)')
title('Signal w-PWM(t)')









