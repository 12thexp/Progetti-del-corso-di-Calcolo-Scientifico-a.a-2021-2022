clearvars, close all
% esegue segment_image a k fissato, disegna l'immagine e il grafico 3D dei
% k colori usati nella quantizzazione

img = imread("mandrill.png"); %immagine originale
imgs = single(img) ./ 255;
imgs = reshape(imgs, [], 3);

%% parametri 
k = 14;
maxiter = 100;
stype = 'plus';

%% segment_image
quanti = segment_image(img, k, maxiter, stype);
quant = single(quanti);
quants = reshape(quant, [], 3);

%% immagine quantizzata
figure
image(quanti)
axis square;
set(gca, 'Visible', 'off')
title("k = " + k + ", start = " + stype)

%% RBG plot dei colori
figure
ax1 = nexttile;
axis(ax1, [0 255 0 255 0 255])
scatter3(quants(:,1),quants(:,2),quants(:,3), 100, quants ./ 255, 'filled')
set(gca,'XLim',[0 255],'YLim',[0 255],'ZLim',[0 255])
xlabel('R', 'Color', 'r')
ylabel('B', 'Color', 'b')
zlabel('G', 'Color', 'g')



