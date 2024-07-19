clearvars, close all;

% soluzione alle DF del pb 
% -D*c" + u*c' + Kc = 1
% u(-L) = 0, u(L) = 0

% calcoliamo due sol approssimate, z1 sull'intervallo (-L, 0), e z2
% sull'intervallo (0, L), per poi unirle in un unico vettore
% questo perche' la sol approssimata su (-L,0) presenta il fenomeno 
% dello strato limite per D -> 0, e va gestito separatamente

% usiamo sempre una discretizzazione di tipo upwind

L = 100;
u = 1; M = 2; K = 1;
Dv = [1 1e-1 1e-2];

n = 200; h = L/n; % passo di griglia
dim = n+2; % nr nodi totale per ciascun intervallo 

f = @(x)1+0.*x;

for i = 1:numel(Dv) % ciclo al variare di D 
    D = Dv(i);
    Pe(i) = (u*h)/(2*D);
    
    Dh = D*(1 + Pe(i)); % discr upwind del termine di trasporto = viscosità artificiale 
    
    % A e' la stessa su entrambi gli intervalli
    a1 = -Dh/h^2 + u/(2*h);
    b1 = 2*Dh/h^2 + K;
    c1 = -Dh/h^2 - u/(2*h);

    d0 = b1*ones(n,1);
    dup = a1*ones(n-1,1);
    dlow = c1*ones(n-1,1);
    A = diag(d0, 0) + diag(dup, 1) + diag(dlow, -1);
    
    
    % ********** INTERVALLO 1 (-L, 0) **********
    xmesh1 = linspace(-L,0,n+2); 
    zsx = 0; zdx = M/(sqrt(u^2+4*D*K));    %BCs
    
    %risolvo il sistema Az1 = v
    f = M*zeros(n,1);
    f(1) = f(1) - c1*zsx;
    f(end) = f(end) - a1*zdx;
    
    z1 = A\f;
    z1 = [zsx; z1; zdx];
    
    % soluzione esatta su (-L, 0)
    l_1 = (u + sqrt(u^2+4*D*K))/(2*D); % > 0
    uex1 = @(x) (M*exp(l_1*x))/(sqrt(u^2+4*D*K));
    
    err1 = norm(uex1(xmesh1)' - z1,'inf');
    %-----------------------------------------------------------
    
    
    % ********** INTERVALLO 2 (0, L) **********
    xmesh2 = linspace(0, L, n+2); 
    zsx = M/(sqrt(u^2 + 4*D*K)); zdx = 0;     %BCs
    
    %risolvo il sistema Az2 = v
    %f = (h^2/Dh)*ones(n,1);
    f = M*zeros(n,1);
    f(1) = f(1) - c1*zsx;
    f(end) = f(end) - a1*zdx;
    
    z2 = A\f;
    z2 = [zsx; z2; zdx];
    
    % soluzione esatta su (0, L)
    l_2 = (-u + sqrt(u^2+4*D*K))/(-2*D); % < 0
    uex2 = @(x) (M*exp(l_2*x))/(sqrt(u^2+4*D*K));
    
    err2 = norm(uex2(xmesh2)' - z2,'inf');
    %-----------------------------------------------------------
    
    
    % ********** SOL COMPLESSIVA **********
    xmesh = [xmesh1 xmesh2];
    uex = [uex1(xmesh1) uex2(xmesh2)];
    z = [z1; z2];
    figure
    plot(xmesh, uex,'b-', xmesh, z,'r-');
    xlim([-10 10])
  
    legend('esatta', 'approssimata')
    err(i,:) = [err1 err2];
end




