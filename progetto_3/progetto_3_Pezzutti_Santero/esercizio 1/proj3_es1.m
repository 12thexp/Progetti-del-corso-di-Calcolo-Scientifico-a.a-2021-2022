clearvars, close all

%PROGETTO 3- ESERCIZIO 1

M = 60;
dx = 1/M;
u0 = zeros(M+1,1);
b = 1; %costante di convezione 
D = 1; %costante di diffusione
xx = linspace(0,1,M+1);
%calcoliamo dt
dtcr = dx^2/2;
dt = 0.8*dtcr;

%calcoliamo C e Pe
C = b*dt/dx;
Pe = b*dx/(2*D);

%soluzione esatta
f = @(t,x)exp(-4*pi^2*1.*t)*sin(2*pi*(x-b.*t));
%imponiamo condizioni iniziali
u0 = f(0,xx);
u0  =  u0';

a1 = (-b/(2*dx)+D/(dx)^2)*dt;
b1 = -2*dt*D/(dx)^2+1;
c1 = (b/(2*dx)+D/(dx)^2)*dt;

uold = u0;
unew = zeros(M+1,1);

kmax  =  2000;
for k = 1:kmax  %ciclo sul tempo 
  tk  =  k*dt;
  
    for j = 2:M %ciclo sullo spazio
        unew(j) = a1*uold(j+1)+b1*uold(j)+c1*uold(j-1);
    end
    
    %imponiamo le condizioni di Dirichlet in x=0,x=1
    unew(1) = f(tk,0);
    unew(end) = f(tk,1);
    
    uold = unew;
end

%soluzione discreta al tempo finale
plot(xx,unew, 'b'), hold on
%soluzione esatta al tempo finale
plot(xx,f(tk,xx), 'r')
title('Soluzione approssimata e soluzione esatta')
legend('soluzione approssimata', 'soluzione esatta')
set(gca, 'Fontsize', 14)
grid on