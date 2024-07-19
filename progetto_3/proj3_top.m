clearvars, close all
%discretizzazione alle differenze finite del BVP
% Ut-lapl(U(x,y)) = f(x,y,t), u = g sul bordo

%soluzione esatta del problema
uex = @(x,y,t) exp(-(y.^2+x.^2)./(1+4.*t))./(sqrt(1+4.*t));
f = @(x,y,t) (2*exp(-(y.^2+x.^2)./(1+4.*t)))./(1+4.*t).^(3/2); % sorgente

a = 0; b = 1;
n = 50; %nodi tot, n-1 intervalli
nint_tot = (n-2)^2; %nr nodi interni totali
nx = n-2;
h = (b-a)/(n-1);

dtcr=h^2/2; %dt critico, ricavato dalla CFL condition
dt=0.8*dtcr;

T = 10;
N = 1000;
%dt = T / N
dt = 1e-8;
mu = dt/(h^2);


%BCs
gS = @(x,y,t)exp(-(x.^2+0*y.^2)./(1+4.*t))/(sqrt(1+4.*t)); %sud
gN = @(x,y,t)exp(-(x.^2+1)./(1+4.*t))/(sqrt(1+4.*t));   %nord
gO = @(x,y,t)exp(-(0*x.^2+y.^2)./(1+4.*t))/(sqrt(1+4.*t)); %ovest
gE = @(x,y,t)exp(-(1+y.^2)./(1+4.*t))/(sqrt(1+4.*t)); %est

z = linspace(0,1,n);  %n punti, n-1 intervalli di ampiezza h  
w = linspace(0,1,n);
[Z,W] = meshgrid(z,w);  %griglia completa

x = z(2:n-1);
y = w(2:n-1);
[X,Y] = meshgrid(x,y);  %griglia punti interni

%matrice A
A = diag((1+2*mu)*ones(1,nx)) + diag(-mu*ones(1,nx-1),-1) + diag(-mu*ones(1,nx-1),1);

u_half = zeros(nx); %matrice valori approssimati a t(k+1/2)
u_new = zeros(nx);  %matrice valori approssimati a t(k+1)

uex0 = @(x,y) uex(x,y,0); %sol esatta a t=0
u_k = uex0(X,Y);    %matrice condizioni iniziali


%ciclo sul tempo
for k = 0:N
    tk = k*dt;  %t(k)
    tk05 = (k+1/2)*dt; %t(k+1/2)
    tk1 = (k+1)*dt; %t(k+1)
    
    %matrice BCs al tempo t(k)
    U = zeros(n);
    U(:,1) = gO(0,W(:,1),tk);
    U(:,end) = gE(1, W(:,1),tk);
    U(1,:) = gS(Z(1,:), 1, tk);
    U(end,:) = gN(Z(1,:), 0, tk);
    
    %matrice BCs al tempo t(k+1/2)
    U_half = zeros(n);
    U_half(:,1) = gO(0,W(:,1),tk05);
    U_half(:,end) = gE(1, W(:,1),tk05);
    U_half(1,:) = gS(Z(1,:), 1, tk05);
    U_half(end,:) = gN(Z(1,:), 0, tk05);
    
    %matrice BCs al tempo t(k)
    U_one = zeros(n);
    U_one(:,1) = gO(0,W(:,1),tk);
    U_one(:,end) = gE(1, W(:,1),tk);
    U_one(1,:) = gS(Z(1,:), 1, tk);
    U_one(end,:) = gN(Z(1,:), 0, tk);
   
    %equazione 1
    b = zeros(nx);
    bcs = zeros(1, nx);
    for j = 1:nx
        for i = 1:nx
            if j == 1
                b(i,j) = (1-2*mu)*u_k(i,j) + mu*u_k(i,j+1) + dt/2*f(Z(i),W(j),tk05);
            elseif j == nx
                b(i,j) = mu*u_k(i,j-1) + (1-2*mu)*u_k(i,j) + dt/2*f(Z(i),W(j),tk05);
            else
                b(i,j) = mu*u_k(i,j-1) + (1-2*mu)*u_k(i,j) + mu*u_k(i,j+1) + dt/2*f(Z(i),W(j),tk05);
            end
        end
        
            if (j == 1)
                bcs = U(2:n-1,1);
                bcs(1) = bcs(1) + U_half(1,2);
                bcs(end) = bcs(end) + U_half(end,2);
            elseif (j == nx)
                bcs = U(2:n-1,end);
                bcs(1) = bcs(1) + U_half(1,end-1);
                bcs(end) = bcs(end) + U_half(end,end-1); 
            else
                bcs = zeros(nx, 1);
                bcs(1) = U_half(1,j);
                bcs(end) = U_half(end,j);
            end
  
         u_half(:,j) = thomas(A,b(:,j) + mu*bcs);
    end
%     figure, surf(X,Y,u_half)
%     title('half')
    
    b = zeros(nx);
    bcs = zeros(1, nx);
    %equazione 2
    for i = 1:nx
        for j = 1:nx
            if i == 1
                b(i,j) = (1-2*mu)*u_half(i,j) + mu*u_half(i+1,j) + dt/2*f(Z(i),W(j),tk05);;
            elseif i == nx
                b(i,j) = mu*u_half(i-1,j) + (1-2*mu)*u_half(i,j) + dt/2*f(Z(i),W(j),tk05);
            else
                b(i,j) = mu*u_half(i-1,j) + (1-2*mu)*u_half(i,j) + mu*u_half(i+1,j) + dt/2*f(Z(i),W(j),tk05);
            end
                  
        end
        
           if (i == 1)
                bcs = U_half(1,2:n-1);
                bcs(1) = bcs(1) + U_one(2,1);
                bcs(end) = bcs(end) + U_one(2,end);
            elseif (i == nx)
                bcs = U_half(end,2:n-1);
                bcs(1) = bcs(1) + U_one(end-1,1);
                bcs(end) = bcs(end) + U_one(end-1,end); 
            else
                bcs = zeros(1, nx);
                bcs(1) = U_one(i,1);
                bcs(end) = U_one(i,end);
           end
        
         u_new(i,:) = thomas(A,b(i,:) + mu*bcs);
         u_k = u_new;
    end
%     figure, surf(X,Y,u_new)
%     title('solo nodi interni')
%     
%     figure, surf(Z, W,uex(Z, W, k));
%     title('sol esatta')
end

%calcolo l'errore
%err = norm(reshape(u_new,[],1) - reshape(uex(X,Y,tk1),[],1))
err_inf = norm(reshape(u_new,[],1) - reshape(uex(X,Y,tk1),[],1),inf)

%creo la matrice completa
U_one(2:n-1, 2:n-1) = u_new; %inserisco i punti interni nella matrice delle Bcs a t(k+1)

surf(Z,W,U_one); 
title('sol approssimata')

figure, surf(Z, W, uex(Z, W, tk1));
title('sol esatta')



