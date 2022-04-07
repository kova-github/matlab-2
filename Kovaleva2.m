%ФАЗОВАЯ МОДУЛЯЦИЯ(ФМ) ВАРИАНТ 1
clc 
clear all
close all 
nfig = 1;
figure('position',[100,100,700,500]);
f0 = 1200;
Vmod = 600;
Vinf = 600;
T0 = 1 / f0;
T = 1 / Vmod;
q = 2 ^ round(Vinf * T);
A = 1;
L0 = 50;
dt = T0 / L0;   
t = 0 : dt : T;
si = zeros(q, length(t));
E = zeros(q, 1);

for i = 0:q-1
    si(i + 1, :) = A*cos(2*pi*f0*t - 2*pi*i/q);
    figure(nfig);
    subplot(2,1,i+1);
    plot(t, si(i + 1,:),'LineWidth',2)
    xlabel('t');
    ylabel('si(t)');
    title(['Сигнал i=', num2str(i)]);
    legend('f0 = 1200,Vmod = 600,Vinf = 600'); 
    axis([0 T -1.5 1.5]);
    grid on
    
end
for i = 0:q-1
    E(i+1) = sum(si(i + 1, :).^2) *dt;
    fprintf('E%d = %d\n', i+1, E(i+1)); 
end

% Задание номер 2. Спектры
nfig = 2;
dF = 1/(10*T);
F = 0:dF:6500;
S = zeros(q, length(F));

for i = 0:q-1
    S(i+1,:) = A*(T/2)*cos(2*pi*i/q)*(sinc((F - f0)*T)+sinc((F + f0)*T)).*exp(-1j*pi*F*T)...
        +(1/1j)*A*(T/2)*sin(2*pi*i/q)*(sinc((F - f0)*T)+sinc((F + f0)*T)).*exp(-1j*pi*F*T);
    figure(nfig); 
    subplot(2,1,i+1);
    plot(F, abs(S(i + 1, :)),'m','LineWidth', 2);
    hold on;
    grid on;
    xlabel('f');
    ylabel('S(f)');
    title(['Спектр сигнала Si(t) дискретной ФМ при i = ', num2str(i) + 1]);
    legend('f0 = 1200,Vmod = 600,Vinf = 600'); 
end

 nfig=3;
I=["1,1,1,1,0","0,0,1,1,0"];
ss = zeros(2, length(F));
ss(2, :) =A*(T/2)*(3*(sinc((F-f0)*T)+sinc((F+f0)*T))-2*(sinc((F-f0)*T)+sinc((F+f0)*T))).*exp(-1j*2*pi*F*T);
ss(1, :) = A*(T/2)*((sinc((F-f0)*T)+sinc((F+f0)*T))-4*(sinc((F-f0)*T)+sinc((F+f0)*T))).*exp(-1j*2*pi*F*T);
   for i=0:q-1
   figure(nfig);
   subplot(2,1,i+1)
    plot(F,abs(ss(i+1, :)),'LineWidth',2);
    xlabel('F');
    ylabel('si(F)');
   title(['спектр последовательностей сигнала при N=5, I=',I(i+1)]);
    legend('f0 = 1200,Vmod = 600,Vinf = 600'); 
    grid on
   end
    
 nfig=4;
ss = zeros(2, length(F));
ss(2, :) = -1*A*(T/2)*((sinc((F-f0)*T)+sinc((F+f0)*T))).*exp(-1j*pi*F*T);
ss(1, :) = A*(T/2)*((sinc((F-f0)*T)+sinc((F+f0)*T)));
   for i=0:q-1
   figure(nfig);
   subplot(2,1,i+1)
    plot(F,abs(ss(i+1, :)),'LineWidth',2);
    xlabel('F');
    ylabel('si(F)');
   title(['Амплитудный спектр отрезка гармоники сигнала при i=', num2str(i)]);
    legend('f0 = 1200,Vmod = 600,Vinf = 600'); 
    grid on
   end
    
    nfig=5;
I=["0,0","1,1","0,1"];
ss = zeros(3, length(F));
ss(2, :) =A*(T/2)*(2*(sinc((F-f0)*T)+sinc((F+f0)*T))).*exp(-1j*2*pi*F*T);
ss(3, :) =A*(T/2)*(-2*(sinc((F-f0)*T)+sinc((F+f0)*T))).*exp(-1j*2*pi*F*T);
ss(1, :) = A*(T/2)*((sinc((F-f0)*T)+sinc((F+f0)*T))-(sinc((F-f0)*T)+sinc((F+f0)*T))).*exp(-1j*2*pi*F*T);
   for i=0:q
   figure(nfig);
   subplot(3,1,i+1)
    plot(F,abs(ss(i+1, :)),'LineWidth',2);
    xlabel('F');
    ylabel('si(F)');
   title(['спектр последовательностей сигнала при N=2, I=',I(i+1)]);
    legend('f0 = 1200,Vmod = 600,Vinf = 600'); 
    grid on
   end

    
    nfig=6;
I=["0,0,1,1,0,0,0,1,1,0","1,1,1,1,1,1,1,1,1,1"];
ss = zeros(2, length(F));
ss(2, :) = A*(T/2)*(-10*(sinc((F-f0)*T)+sinc((F+f0)*T))).*exp(-1j*2*pi*F*T);
ss(1, :) = A*(T/2)*(6*(sinc((F-f0)*T)+sinc((F+f0)*T))-4*(sinc((F-f0)*T)+sinc((F+f0)*T))).*exp(-1j*2*pi*F*T);
   for i=0:q-1
   figure(nfig);
   subplot(2,1,i+1)
    plot(F,abs(ss(i+1, :)),'LineWidth',2);
    xlabel('F');
    ylabel('si(F)');
   title(['спектр последовательностей сигнала при N=10, I=',I(i+1)]);
    legend('f0 = 1200,Vmod = 600,Vinf = 600'); 
    grid on
   end
