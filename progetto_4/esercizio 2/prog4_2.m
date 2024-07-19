clearvars, close all

%PROGETTO 4- ESERCIZIO 2

%definisco la matrice K
K=[1 0; 0 1/(2^10)];
%definisco il vettore y
y=[1;2^(-10)];

%calcolo soluzione esatta di Kx=y
x=K\y;

%problema 1 perturbato
p1=[0;2^(-10)];
yp1=y+p1;
%soluzione perturbata
xp1=K\yp1;
r1=norm(xp1-x)/norm(p1);

%problema 2 perturbato
p2=[2^(-10);0];
yp2=y+p2;
%soluzione perturbata
xp2=K\yp2;
r2=norm(xp2-x)/norm(p2);

%calcolo il condizionamento della matrice K
condK=cond(K); 

alpha=[0.001:0.001:1-0.001];

%problema perturbato 1
for i=1:numel(alpha)
    a=alpha(i);
    %definisco la matrice K in funzione di a
    Ka=[1 0; 0 (1/(1-a)^10)*(1/2^(10))];
    %calcolo la soluzione del sistema Ka*xp1a=yp1
    xp1a=Ka\yp1;
    %calcolo la quantità r1 in funzione di a
    r1a(i)=norm(xp1a-x)/norm(p1);

end

%trovo il valore di alpha che minimizza r1a
[r1min,i1min] = min(r1a);
alpha1min=alpha(i1min);

%calcolo condizionamento della matrice avente come alpha=alpha1min 
%e la soluzione perturbata
Ka1=[1 0; 0 (1/(1-alpha1min)^10)*(1/2^(10))];
cond1=cond(Ka1);
xp1a1=Ka1\yp1;



%problema perturbato 2
for i=1:numel(alpha)
    a=alpha(i);
    %definisco la matrice K in funzione di a
    Ka=[1 0; 0 (1/(1-a)^(10))*(1/2^(10))];
    %calcolo la soluzione del sistema Ka*xp2a=yp2
    xp2a=Ka\yp2;
    %calcolo la quantità r2 in funzione di a
    r2a(i)=norm(xp2a-x)/norm(p2);

end
%trovo il valore di alpha che minimizza r2a
[r2min,i2min] = min(r2a);
alpha2min=alpha(i2min);

%calcolo condizionamento della matrice avente come alpha=alpha2min 
%e la soluzione perturbata
Ka2=[1 0; 0 (1/(1-alpha2min)^10)*(1/2^(10))];
cond2=cond(Ka2);
xp2a2=Ka2\yp2;



%calcolo il condizionamento di Ka
for i=1:numel(alpha)
    a=alpha(i);
    Ka=[1 0; 0 (1/(1-a)^10)*(1/2^(10))];
    %calcolo il condizionamento della matrice Ka
    condKa(i)=cond(Ka);
end

%trovo il valore di alpha che minimizza condKa
[condKamin,icondKamin]=min(condKa);
alphacmin=alpha(icondKamin);


