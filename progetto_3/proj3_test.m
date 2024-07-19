clearvars, close all
%discretizzazione alle differenze finite del BVP
% Ut-lapl(U(x,y)) = f(x,y,t), u = g sul bordo

%soluzione esatta del problema calcolata con separazione variabili
uex = @(x,y,t) exp(-(y.^2+x.^2)./(1+4.*t))./(sqrt(1+4.*t));
sorgente = @(x,y,t) (2*exp(-(y.^2+x.^2)./(1+4.*t)))./(1+4.*t).^(3/2); % sorgente

a = 0; b = 1;
n = 10;  %nr intervalli
nint = (n-1)^2; %nr nodi interni 
h = (b-a)/n;
dt = 1e-2;

mu = dt/(h^2);


%BCs
gN = @(x,y,t)exp(-(x.^2+0*y.^2)./(1+4.*t))/(sqrt(1+4.*t));
gS = @(x,y,t)exp(-(x.^2+1)./(1+4.*t))/(sqrt(1+4.*t));
gO = @(x,y,t)exp(-(0*x.^2+y.^2)./(1+4.*t))/(sqrt(1+4.*t)); 
gE = @(x,y,t)exp(-(1+y.^2)./(1+4.*t))/(sqrt(1+4.*t)); 

z = linspace(0,1,n+1);  %n+1 punti, n intervalli di ampiezza h  
w = linspace(0,1,n+1);
[Z,W] = meshgrid(z,w);  %griglia completa

x = z(2:n);
y = w(2:n);
[X,Y] = meshgrid(x,y);  %griglia punti interni

%matrice u_k delle condizioni iniziali
uex0 = @(x,y) uex(x,y,0); %sol esatta a t=0
u_kv = reshape(uex0(X,Y)',[],1); % cond iniz, scorro le righe

%matrice A
blockA = diag((1+mu)*ones(1,n-1)) + diag(-mu/2*ones(1,n-2),-1) + diag(-mu/2*ones(1,n-2),1);
Ac = repmat({blockA}, 1, n-1);  
A = blkdiag(Ac{:}); %creo A a blocchi

%A = diag((1+2*mu)*ones(1,nint)) + diag(-mu*ones(1,nint-1),-1) + diag(-mu*ones(1,nint-1),1);

%matrice B
B = diag((1-2*mu)*ones(1,nint)) + diag(mu*ones(1,nint-n),-(nint-2)) + diag(mu*ones(1,nint-n),nint-2);



for k = 1:10
    
    %calcolo f(x,y,k+1/2)
    f = @(x,y) sorgente(x,y,k+1/2); %f al tempo k+1/2
    fv = reshape(f(X,Y),[],1); %valuto f sui punti interni, scorro f per colonne
    
    
    %% BCs 2
    north = gN(Z(1,2:n),1,k+1/2);
    north(1) = north(1) + gE(Z(2),W(1),k+1/2);
    north(end) = north(end) + gE(Z(2),W(n+1),k+1/2);
    
    if n+1 > 4
        for i = 1:n-3
            r(:,i) = [gO(Z(1), W(i+1), k+1/2); zeros(n-3,1); gE(Z(n+1), W(i+1), k+1/2)];
        end
    else
        r = 0;
    end
    
    south = gS(Z(1,2:n),0,k+1/2);
    south(1) = south(1) + gE(Z(n),W(1),k+1/2);
    south(end) = south(end) + gE(Z(n),W(n+1),k+1/2);
    
    bcs2 = [north'; reshape(r,[],1); south'];
    
    %% BCs 1
    east = gE(1,W(1,2:n),k+1);
    east(1) = east(1) + gN(Z(n),W(1),k+1);
    east(end) = east(end) + gS(Z(n),W(n),k+1);
    
    if n+1 > 4
        for i = 1:n-3
            r(:,i) = [gN(Z(i+1), W(1), k+1); zeros(n-3,1); gS(Z(i+1), W(n+1), k+1)];
        end
    else
        r = 0;
    end
    
    west = gO(0,W(n,2:n),k+1);
    west(1) = west(1) + gN(Z(2),W(1),k+1);
    west(end) = west(end) + gS(Z(2),W(n+1),k+1);
    
    
    bcs1 = [west'; reshape(r,[],1); east'];
    
   
    %%
    b = B*u_kv + mu*bcs1 + 1/2*dt*fv;    %eq 1
    u_half = thomas(A,b);
%     u_k = reshape(u_kv, [], sqrt(numel(u_kv)));
%     for j = 2:n-1
%         b(j) = u_k(j,:) + mu*(u_k(j-1,:) - 2*u_k(j,:) + u_k(j+1,:)) + dt*1/2*f(X(j),Y(j),k+1/2);
%         
%         u_new(j,:)=(thomas(A,b))'; 
%     end
    
    b = B*u_half + mu*bcs2 + 1/2*dt*fv;  %eq 2
    u_new = thomas(A,b);
    
    u_kv = u_new;

end


u_full = zeros(n+1);  %creo una matrice per la sol approssimata su tutta la griglia
l = 1;
for i = 2:n
    for j = 2: n
        u_full(i,j) = u_kv(l);  %inserisco i punti interni
        l=l+1;
    end
end

% u_full(:,1) = gN(z,1,k+1);  %inserisco i bordi
% u_full(:,n+1) = gS(z,0,k+1);
% u_full(1,:) = gO(0,w,k+1)';
% u_full(n+1,:) = gE(1,w,k+1)';

u_full(:,1) = gN(z,1,0);  %inserisco i bordi
u_full(:,n+1) = gS(z,0,0);  %qui prima gS e gN erano invertite, ma forse dovevo cambiare i nomi delle funz originali
u_full(1,:) = gO(0,w,0)';
u_full(n+1,:) = gE(1,w,0)';

Uex = uex(Z,W,k); %calcolo la sol esatta sulla griglia

err = norm(reshape(u_full,[],1) - reshape(Uex,[],1),inf)

surf(Z,W,u_full); 
title('sol approssimata con condizioni al bordo')

figure, surf(Z, W,uex(Z, W, 0));
title('sol esatta')

u_k = reshape(u_kv, n-1, n-1);
figure, surf(X,Y,u_k)
title('solo nodi interni')



