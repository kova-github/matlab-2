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
dt = T0 / L0; 
df=f0/n;
t = 0 : dt : T;
s = zeros(q, length(t));
L0 = 50;
 

for i = 0:q-1
    Sn = A*cos(2*pi*f0*t - 2*pi*i/q);
    s(i+1,:) = Sn;
end

baz1 = sqrt(2/T)*cos(2*pi*f0*t);
baz1 = baz1/norm(baz1);
    
baz2 = sqrt(2/T)*sin(2*pi*f0*t);
baz2 = baz2/norm(baz2);    

fprintf('произведение базисных функций f(1) и f(1) = %1.f \n',round(baz1*baz1));
fprintf('произведение базисных функций f(1) и f(2) = %1.f \n',round(baz1*baz2));
fprintf('произведение базисных функций f(2) и f(2) = %1.f \n',round(baz2*baz2));

s1 = zeros(1, length(s(i+1,:)));
s2 = zeros(1, length(s(i+1,:)));  

for i = 0:q
    s1 = s(i,:)*baz1';
    s2 = s(i,:)*baz2';
    
    plot(s1,s2,'k');
    holf on;
    grid on;
    
end

x = -6:6;
y = x;
plot(x,y);
hold on;
grid on;
axis('square');