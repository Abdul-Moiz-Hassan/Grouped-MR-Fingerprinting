clc
clear all
close all

T2DisplayMax=2000;
T1DisplayMax=3000;

load T1_CPU.txt
load T2_CPU.txt

imT1 = reshape(T1_CPU,[240,240]);
imT2 = reshape(T2_CPU,[240,240]);

load t1big.mat;
load t2big.mat;

subplot(231); imagesc(imT1,[0 T1DisplayMax]); colormap(jet); axis image; axis off; title('GPU T_1'); colorbar
subplot(234); imagesc(imT2,[0 T2DisplayMax]); colormap(jet); axis image; axis off; title('GPU T_2'); colorbar

subplot(232); imagesc(t1big,[0 T1DisplayMax]); colormap(jet); axis image; axis off; title('Actual T_1'); colorbar
subplot(235); imagesc(t2big,[0 T2DisplayMax]); colormap(jet); axis image; axis off; title('Actual T_2'); colorbar

d1=abs(t1big-imT1);
d2=abs(t2big-imT2);

subplot(233); imagesc(d1,[0 T1DisplayMax]); colormap(jet); axis image; axis off; title('Difference T_1'); colorbar
subplot(236); imagesc(d2,[0 T2DisplayMax]); colormap(jet); axis image; axis off; title('Difference T_2'); colorbar

snr1 = calc_snr(imT1,t1big)
snr2 = calc_snr(imT2,t2big)
