clear all
close all
clc

%% Seleção de antenas (AS), MRC e Alamouti.

Nbits = 200e3;
SNR_M = 20;
fd = 100; %freq Doppler
M = 2;
tb = 1/Nbits;

info = randi([0 1], 1, Nbits);
info_mod = pskmod(info, M);

%% AS (1Tx, 2Rx)

canal_as1 = rayleighchan(tb, fd);
canal_as1.StoreHistory = 1;
sinal_as1 = filter(canal_as1, info_mod);
ganho_canal_as1 = canal_as1.PathGains;

canal_as2 = rayleighchan(tb, fd);
canal_as2.StoreHistory = 1;
sinal_as2 = filter(canal_as2, info_mod);
ganho_canal_as2 = canal_as2.PathGains;

%% MRC (1Tx, 2Rx)

canal_mrc1 = rayleighchan(tb, fd);
canal_mrc1.StoreHistory = 1;
sinal_mrc1 = filter(canal_mrc1, info_mod);
ganho_canal_mrc1 = canal_mrc1.PathGains;

canal_mrc2 = rayleighchan(tb, fd);
canal_mrc2.StoreHistory = 1;
sinal_mrc2 = filter(canal_mrc2, info_mod);
ganho_canal_mrc2 = canal_mrc2.PathGains;


%% Alamouti (2Tx, 1Rx)

tx1 = zeros(1, Nbits);
tx2 = zeros(1, Nbits);

tx1(1:2:end) = info_mod(1:2:end);
tx2(1:2:end) = info_mod(2:2:end);
info_mod_tx1 = info_mod(2:2:end);
info_mod_tx2 = info_mod(1:2:end);
tx1(2:2:end) = -conj(info_mod_tx1);
tx2(2:2:end) = conj(info_mod_tx2);

canal_al1 = rayleighchan(tb, fd);
canal_al1.StoreHistory = 1;
sinal_al1 = filter(canal_al1, tx1);
ganho_canal_al1 = transpose(canal_al1.PathGains);

canal_al2 = rayleighchan(tb, fd);
canal_al2.StoreHistory = 1;
sinal_al2 = filter(canal_al2, tx2);
ganho_canal_al2 = transpose(canal_al2.PathGains);

for SNR = 0:SNR_M
    
    % técnica de seleção de antenas
    rx_as1 = awgn(sinal_as1, SNR, 'measured');
    rx_as1_eq = rx_as1.*ganho_canal_as1';
    rx_as2 = awgn(sinal_as2, SNR, 'measured');
    rx_as2_eq = rx_as2.*ganho_canal_as2';
    rx_as = max(rx_as1_eq, rx_as2_eq);   % escolhe o melhor sinal recebido fazendo o maximo entre os dois
    
    % técnica de multiplas antenas
    rx_mrc1 = awgn(sinal_mrc1, SNR, 'measured');
    rx_mrc1_eq = rx_mrc1.*ganho_canal_mrc1';
    rx_mrc2 = awgn(sinal_mrc2, SNR, 'measured');
    rx_mrc2_eq = rx_mrc2.*ganho_canal_mrc2';
    rx_mrc = rx_mrc1_eq + rx_mrc2_eq;    % Soma os dois sinais
      
    % técnica alamout 
    sinal_al = sinal_al1 + sinal_al2;
    rx_alamouti = awgn(sinal_al, SNR, 'measured');
    al_S = zeros(1, length(info));
    al_S(1:2:end) = conj(ganho_canal_al1(1:2:end)).*rx_alamouti(1:2:end) + ganho_canal_al2(2:2:end).*conj(rx_alamouti(2:2:end));
    al_S(2:2:end) = conj(ganho_canal_al2(1:2:end)).*rx_alamouti(1:2:end) - ganho_canal_al1(2:2:end).*conj(rx_alamouti(2:2:end));

    info_demod_as = pskdemod(rx_as, M);
    info_demod_mrc = pskdemod(rx_mrc, M);
    info_demod_al = pskdemod(al_S, M);
    
    [err_as(SNR + 1), taxa_as(SNR + 1)] = biterr(info, info_demod_as);
    [err_mrc(SNR + 1), taxa_mrc(SNR + 1)] = biterr(info, info_demod_mrc);
    [err_al(SNR + 1), taxa_al(SNR + 1)] = biterr(info, info_demod_al);
end

figure(1)
semilogy([0:SNR_M], taxa_as); hold on;
semilogy([0:SNR_M], taxa_mrc);
semilogy([0:SNR_M], taxa_al);
title('AS(1Tx, 2Rx) x MRC(1Tx, 2Rx) x Alamouti(2Tx, 1Rx)');
legend('AS', 'MRC', 'Alamouti');
