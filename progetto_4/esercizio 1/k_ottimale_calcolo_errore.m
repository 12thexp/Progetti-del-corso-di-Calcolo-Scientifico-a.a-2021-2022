clearvars, close all
% calcola l'errore per k = 1:kmax, e fa il grafico dell'andamento
% dell'errore

img = imread("mandrill.png");
imgs = single(img) ./ 255;
imgs = reshape(imgs, [], 3);
kmax = 8;
maxiter = 100;
stype = 'plus';

%% calcolo k ottimale
for k = 1:kmax
    
    quanti = segment_image(img, k, maxiter, stype);
    quant = single(quanti);
    quants = reshape(quant./255, [], 3);

    %% errore quadratico medio
    err = 0;
    n = 512^2;
    for i = 1:n
        d = norm(imgs(i,:)) - norm(quants(i,:)); %differenza tra i valori del colore orig. e approx. per l'iesimo pixel
        cell_error = 1/n*(d^2);
        err = err + cell_error;
    end
    ERR(k) = err;
end

%% grafico dell'errore
k_v = 1:1:k;
plot(k_v, ERR, 'o-')
xlabel('k')
ylabel('error')
