clear all
close all
clc


Nbits = 200e3;
SNR_M = 20;
fd = 100; %freq Doppler
M = 2;
tb = 1/Nbits;

info = randi([0 1], 1, Nbits);
info_mod = pskmod(info, M);

%% SISO (1 Tx, 1 Rx)

canal_siso = rayleighchan(tb, fd);
canal_siso.StoreHistory = 1;
sinal_siso= filter(canal_siso, info_mod);
ganho_canal_siso = canal_siso.PathGains;



%% Alamouti (2Tx)

tx1 = zeros(1, Nbits);
tx2 = zeros(1, Nbits);

tx1(1:2:end) = info_mod(1:2:end);
tx2(1:2:end) = info_mod(2:2:end);
info_mod_tx1 = info_mod(2:2:end);
info_mod_tx2 = info_mod(1:2:end);
tx1(2:2:end) = -conj(info_mod_tx1);
tx2(2:2:end) = conj(info_mod_tx2);

%% Alamouti 1Rx

canal_al1 = rayleighchan(tb, fd);
canal_al1.StoreHistory = 1;
sinal_al1 = filter(canal_al1, tx1);
ganho_canal_al1 = transpose(canal_al1.PathGains);

canal_al2 = rayleighchan(tb, fd);
canal_al2.StoreHistory = 1;
sinal_al2 = filter(canal_al2, tx2);
ganho_canal_al2 = transpose(canal_al2.PathGains);


%% Alamouti 2Rx

canal_al3 = rayleighchan(tb, fd);
canal_al3.StoreHistory = 1;
sinal_al3 = filter(canal_al3, tx1);
ganho_canal_al3 = transpose(canal_al3.PathGains);

canal_al4 = rayleighchan(tb, fd);
canal_al4.StoreHistory = 1;
sinal_al4 = filter(canal_al4, tx2);
ganho_canal_al4 = transpose(canal_al4.PathGains);

for SNR = 0:SNR_M
    
    % técnica SISO
    rx_siso = awgn(sinal_siso, SNR, 'measured');
    rx_siso_eq = rx_siso.*ganho_canal_siso';
    
    % tecnica Alamouti (2Tx, 1Rx)
    sinal_al = sinal_al1 + sinal_al2;
    rx_alamouti_1 = awgn(sinal_al, SNR, 'measured');
    al_S = zeros(1, length(info));
    al_S(1:2:end) = conj(ganho_canal_al1(1:2:end)).*rx_alamouti_1(1:2:end) + ganho_canal_al2(2:2:end).*conj(rx_alamouti_1(2:2:end));
    al_S(2:2:end) = conj(ganho_canal_al2(1:2:end)).*rx_alamouti_1(1:2:end) - ganho_canal_al1(2:2:end).*conj(rx_alamouti_1(2:2:end));

    
    % tecnica Alamouti (2Tx, 2Rx)
    sinal_a2 = sinal_al3 + sinal_al4;
    rx_alamouti_2 = awgn(sinal_a2, SNR, 'measured');
    al_Sb = zeros(1, length(info));
    al_Sb(1:2:end) = conj(ganho_canal_al1(1:2:end)).*rx_alamouti_1(1:2:end) + ganho_canal_al2(2:2:end).*conj(rx_alamouti_1(2:2:end)) + conj(ganho_canal_al3(1:2:end)).*rx_alamouti_2(1:2:end) + ganho_canal_al4(2:2:end).*conj(rx_alamouti_2(2:2:end));
    al_Sb(2:2:end) = conj(ganho_canal_al2(1:2:end)).*rx_alamouti_1(1:2:end) - ganho_canal_al1(2:2:end).*conj(rx_alamouti_1(2:2:end)) + conj(ganho_canal_al4(1:2:end)).*rx_alamouti_2(1:2:end) - ganho_canal_al3(2:2:end).*conj(rx_alamouti_2(2:2:end));
    
    
    info_demod_siso = pskdemod(rx_siso_eq, M);
    info_demod_al1 = pskdemod(al_S, M);
    info_demod_al2 = pskdemod(al_Sb, M);
    
    [err_siso(SNR + 1), taxa_siso(SNR + 1)] = biterr(info, info_demod_siso);
    [err_al1(SNR + 1), taxa_al1(SNR + 1)] = biterr(info, info_demod_al1);
    [err_al2(SNR + 1), taxa_al2(SNR + 1)] = biterr(info, info_demod_al2);

end

figure(1)
semilogy([0:SNR_M], taxa_siso); hold on;
semilogy([0:SNR_M], taxa_al1);
semilogy([0:SNR_M], taxa_al2);
title('SISO(1Tx, 1Rx) x  Alamouti(2Tx, 1Rx) x Alamouti(2Tx, 2Rx)');
legend('SISO', 'Alamouti 1Rx', 'Alamouti 2Rx');