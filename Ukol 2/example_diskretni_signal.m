%--------------------------------------------------------------------------
%Example Cviceni 2 PMZS - Diskretni signal
%--------------------------------------------------------------------------
clear all;
close all;

%--------------------------------------------------------------------------
%Signal definition w[n]
N=10;            %signal period 
POINT_NUMBER = N; %period numbers
Ts=1/50; %sampling frequency
% Harmonic signal defintion
Wmax = 4;
OMEGA = 10*pi;
phase = pi;

n=0:1:POINT_NUMBER-1;				%time axe
w=Wmax*cos(OMEGA*n*Ts+phase);		% signal w[n]

%Signal graphs w[n]
subplot(2,1,1)
stem(n,w); 
xlabel('n')
ylabel('w[n]')
title('Harmonic signal w[n]')

%Signal graphs w[n]
t1 = Ts.*n;
subplot(2,1,2)
stem(t1,w); 
xlabel('t[s]')
ylabel('w[t]')
title(' Harmonic signal w[t]')

