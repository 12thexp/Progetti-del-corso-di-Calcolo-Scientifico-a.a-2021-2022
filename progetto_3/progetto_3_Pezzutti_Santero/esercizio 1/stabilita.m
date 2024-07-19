clearvars, close all

%%%PROGETTO 3- ESERCIZIO 1 punto 2
%%stabilità del modello proposto

%Pe=x
%C=y

x=linspace(0,20,60);
%Scelgo un valore per k e deltax
k=500;
deltax=1/60;

y=(2.*x-2.*x.*cos(k*deltax))./((1-cos(k*deltax)).^2+(x*sin(k*deltax)).^2);

area(x,y), grid, ylim([-1 3])
xlabel('Pe'), ylabel('C')
title('Metodo stabile')
set(gca, 'Fontsize', 14)