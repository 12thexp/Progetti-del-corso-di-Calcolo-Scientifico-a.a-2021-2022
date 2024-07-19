clearvars, close all

%PROGETTO 3- ESERCIZIO 1

M = 60;
dx = 1/M;
u0 = zeros(M+1,1);
n = 0;
b = 1;
D = 1;
%dt = 0.0001;
xx = linspace(0,1,M+1);
dtcr = dx^2/2; % condizione CFL per schema esplicito 
dt = 0.8*dtcr;

C = b*dt/dx;
Pe = b*dx/(2*D);


f = @(t,x)exp(-4*pi^2*1.*t)*sin(2*pi*(x-b.*t));
u0 = f(0,xx);
u0  =  u0';

 a1 = (-b/(2*dx)+D/(dx)^2)*dt;
 b1 = -2*dt*D/(dx)^2+1;
 c1 = (b/(2*dx)+D/(dx)^2)*dt;

uold = u0;
unew = zeros(M+1,1);

kmax  =  2000;
for k = 1:kmax
  tk  =  k*dt;
  
    for j = 2:M
        unew(j) = a1*uold(j+1)+b1*uold(j)+c1*uold(j-1);
    end
    
    unew(1) = f(tk,0);
    unew(end) = f(tk,1);
    uold = unew;
end

   plot(xx,unew, 'b'), hold on
   plot(xx,f(tk,xx), 'r')