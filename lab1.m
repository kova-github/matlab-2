clc
clear all
nfig = 1;
close all

f0 = 1200; 
Vm = 600; 
Vi = 600; 
T0 = 1 / f0; 
T = 1 / Vm; 
q = 2 ^ round(Vi * T); 

A = 1; 
ot = 64;
dt = T0 / ot;   
t = 0 : dt : T; 
si = zeros(q, length(t));
res = zeros(q,q);

for i = 0:q-1
    si(i + 1, :) = A*cos(2*pi*f0*t - 2*pi*i/q);
 
    figure(nfig); nfig = nfig +1;
    plot(t, si(i + 1, :),'b','LineWidth', 2);
    grid on
    xlabel('t');
    ylabel('s(t)');
    title(['Сигнал ', num2str(i) + 1]);
end


eng1 = sum(si(1,:) .^2)*dt
eng2 = sum(si(2,:) .^2)*dt


%2

dF = 1/(10*T);

F = 0:dF:6500;
S = zeros(q, length(F));

for i = 0:q-1
    S(i+1,:) = A*cos(2*pi*i/q)*(sinc((F - f0)*T)+sinc((F + f0)*T)).*exp(1j*pi*F*T) ...
        + (1/1j)*A*sin(2*pi*i/q)*(sinc((F - f0)*T)+sinc((F + f0)*T)).*exp(1j*pi*F*T);
    figure(nfig); nfig = nfig +1;
    plot(F, abs(S(i + 1, :)),'b','LineWidth', 2);
    hold on;
    grid on;
    xlabel('f');
    ylabel('S(f)');
    title(['Спектр сигнала ', num2str(i) + 1]);
end



N = 5;
multI = [0,0,0,0,0]; 
Si = zeros(1, length(F));
for k = 1:N
    Si = Si + S(multI(k) + 1, :).*exp(-(1j*2*pi*(k-1)*T)*F);
end

figure(nfig); nfig = nfig +1;
plot(F, abs(Si),'m');
hold on;
grid on;
xlabel('f');
ylabel('S(f)');
title('Спектр последовательности');