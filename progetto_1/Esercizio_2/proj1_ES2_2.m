clearvars, close all;

% soluzione alle DF del pb 
% -D*c" + u*c' = 1
% u(-L) = u0, u(L) = 1

L = 1;
x0 = -L; xL = L; % m 
u = 1; % m/s 
Dv = [1 1e-1 1e-2]; % m^2/s

n = 70; h = (xL-x0)/n; % passo di griglia
dim = n+1; % nr nodi totale
xh = linspace(x0, xL, n+1); 

u0 = 0; uL = 1;
f = @(x)1+0.*x;

% soluzione esatta 
uex  =  @(x, D) (1/(exp(2*u/D)+1))*(exp((u/D)*(1+x))+x);

for i = 1:numel(Dv) % ciclo al variare di D 
    D = Dv(i);
    Pe(i) = (u*h)/(2*D);
    
    %Dh = D;  % DF centrate "pure"
    Dh = D*(1+Pe(i)); % discr upwind del termine di trasporto   =  viscosità artificiale 
    
    A = zeros(dim-2);
    f = zeros(dim-2,1);
     
    a1 = -Dh/h^2 + u/(2*h);
    b1 = 2*Dh/h^2;
    c1 = -Dh/h^2 - u/(2*h);

    d0 = b1*ones(dim-2,1);
    dp1 = a1*ones(dim-3,1);
    dm1 = c1*ones(dim-3,1);
    A = diag(d0, 0) + diag(dp1, 1) + diag(dm1, -1);
    f(1) = f(1) - c1*u0;
    f(end) = f(end) - a1*uL;
    
    uh = A\f;
    uh = [u0;uh;uL];
    
    figure
    plot(xh, uex(xh,D), 'b-',xh, uh, 'r-');
    legend('esatta', 'approssimata')
    str_title = strcat('Pe = ',sprintf('%1.0g',Pe(i)));
    title(str_title)
    
    err(i) = norm(uex(xh,D)' - uh,'inf');
    
end
