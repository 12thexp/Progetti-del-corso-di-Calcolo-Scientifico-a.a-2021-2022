clearvars, close all
format short e

%n nodi interni, n+2 nodi totali
nn = [20 40 80 160];
u0 = 1; uL = 1;

%come spiegato nella relazione, usiamo il cambio di variabili u = z+u0; z = u-u0
%ipotizzo che per qualunque u0, uL valga u0 = uL
%di fatto il cambio e' visibile solo sui termini noti
for i = 1:numel(nn)
    n = nn(i);
    h = 2/(n+2-1);
    z = zeros(n,1);
    
    %sol esatta
    xmesh = linspace(-1,1,n+2);
    solinit = bvpinit(xmesh,@guess);
    options = bvpset('RelTol',1e-3, 'AbsTol',1e-6);
    sol = bvp4c(@bvpfcn, @bcfcn, solinit, options);
    
    sol_y = deval(sol, xmesh);  %sol.x non corrisponde ai punti di mesh, quindi ricalcolo la sol sui punti della griglia
    sol_y = sol_y(1, :)'; %traspongo
    
    %La sol. approssimata deve risolvere il sistema: Az - sqrt(Bz)^2 - b = 0;
    A = diag(-2*ones(1,n)) + diag(ones(1,n-1),1) + diag(ones(1,n-1),-1);
    A = A/(h^2);
    
    B = diag(ones(1,n-1),1) - diag(ones(1,n-1),-1);
    B = B/(2*h);
    
    b = zeros(n,1);
    %b(1) = 0/(h^2);
    %b(n) = 0/(h^2);
    
    toll = 1e-4;
    err = 1;
    it = 0;
    itmax = 200;
    
    %Newton
    while (err > toll) && (it < itmax)
        Bz = B*z;
        Bz2 = (Bz).^2;
        Az = A*z;
        f = A*z - sqrt(1+(Bz).^2) - b;
        
        %costruiamo lo Jacobiano di f
        J = zeros(n,n);
        D = zeros(n,n);
        D = diag((h/2)*(Bz)./(sqrt(1+(Bz2))));
        M = D*B;
        J = A+M;
        
        diff = J\-f;
        z_new = diff+z;
        z = z_new;
        err = norm(diff);  %calcoliamo l'errore a ciascun passo
                           %di Newton come distanza dalla sol esatta
        it = it+1;
    end
    
    u_new = z_new + u0;  %torniamo alla variabile u
    u_new = [u0;u_new;uL];  %orliamo con le condizioni al bordo
    
    %grafico per n = 20
    if n == nn(1)
        figure
        plot(xmesh, sol_y,'b*') %sol esatta
        hold on
        plot(xmesh, u_new, 'r-');  %sol approssimata
        grid on
        legend('sol. esatta','sol. approssimata')
        title("n = "+ n)
    end
    
    error(i) = norm(u_new-sol_y,'inf');  %errore finale al variare di n
end

error

for i = 1 : numel(nn)-1
    p = log(error(i)/error(i+1)) / log(2)   %ordine di convergenza
end

%grafico dell'errore
figure
hh = 2./(nn - 1);
loglog(nn, hh.^2, 'black-')    %grafico n vs errore, h^2
grid on
hold on
loglog(nn, error, 'r*', nn, error(end)*(hh./hh(end)).^2, 'b-')  %h^2 traslato e sovrapposto all'errore
set(gca,'FontSize',11)
legend('h^2', 'errore', "h^2 traslato")
xlabel('n')
