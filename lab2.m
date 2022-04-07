%ФАЗОВАЯ МОДУЛЯЦИЯ(ФМ) ВАРИАНТ 1
clc 
clear all
close all 
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
end


% Задание номер 2. Спектры
nfig = 1;
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
    ylabel('Si(f)');
    title(['Спектр сигнала Si(f) дискретной ФМ при i = ', num2str(i)]);
    legend('f0 = 1200,Vmod = 600,Vinf = 600'); 
end

    
 nfig=2;
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

   
nfig=3;
N=1000;
l=round((q-1)*rand(1,N));
ss=zeros(1,length(F));
for k=1:N
    ss=ss+S(l(k)+1,:).*exp(-1j*2*pi*(k-1)*F*T);
end
figure(nfig);
    plot(F,abs(ss),'LineWidth',2);
    xlabel('F');
    ylabel('si(F)');
   title('N=1000');
    legend('f0 = 1200,Vmod = 600,Vinf = 600'); 
    grid on
    

   
   
   