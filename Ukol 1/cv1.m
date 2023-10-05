%% Vykresleni signalu

t = (0:0.01:3);
w = 5*cos(30*t+pi()/4);

plot(t,w);
xlabel('Cas (t)');
ylabel('Amplituda');