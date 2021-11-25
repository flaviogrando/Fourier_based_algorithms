%%%%%%%%%%%%%%%%%%%     TESTES ESTIMADOR DFT    %%%%%%%%%%%%%%%%%%%%%%
% Flavio Lori Grando
% GPGEI - 2019

clc
clear all
close all


set(0, 'defaultAxesFontSize',12);
set(0, 'defaultAxesFontName','times');

% ---------------------------------------------------------------------
% INICIALIZAÇÕES

N = 250;      % número de amostras
m = 250;      % passo de janelamento
Fs = 12500;   % Taxa de amostragem
f0 = 50;      % freq. fundamental (teórica)
noise = 60;    % nível de ruído em dB ( 0 = sem ruído)
T = 100.02;        % tempo total (em segundos)
param = 0;    % freq. de modulação, por exemplo
type = 0;     % seleciona tipo de teste (ver gerador_sinais.m)

freqs = [50.00];        % Vetor de frequências

% phases = [0 0;           % vetor de fases
%          -2*pi/3 0;
%          +2*pi/3 0];    
phases(1,1) = deg2rad(0);
phases(2,1) = deg2rad(-120);
phases(3,1) = deg2rad(+120);

amps = [1 0.3;           % vetor de amplitudes
        1 0.3;
        1 0.3]; 

% ------------------------------------------------------------------------------------------------------------------
% GERAÇÃO DO SINAL
[Va, Vb, Vc, t, refs, mod] = gerador_sinais(f0, freqs, phases, amps,Fs, T,noise, type, param);

% figure, hold on, grid on
% %plot(t,mod*10)
% plot(t,Va)
% plot(t,Vb)
% plot(t,Vc)
% legend('Va','Vb','Vc'), title('Sinais')
% xlabel('tempo (s)'), ylabel('amplitude')

%-------------------------------------------------------------------------------------------------------------------
% SEGMENTAÇÃO (JANELAMENTO)
[Va_seg, Vb_seg, Vc_seg, t_seg, ref_seg] = segmentador(Va, Vb, Vc, t, refs, N, m);

% % Plot amostragem (ciclo a ciclo)
% w = 5;   % numero de janelas no plot
% figure
% for i=1:w
%     subplot(w,1,i), hold on, grid on
%     plot(t_seg(i,:),Va_seg(i,:),'o-')
% end
% 
% figure %Plot referências
% subplot(3,1,1), hold on, grid on
% plot(t_seg(:,1), ref_seg(1,:), 'o-')
% xlabel('tempo (s)'), ylabel('freq (Hz)'), title('Referências (segmentada)')
% subplot(3,1,2), hold on, grid on
% plot(t_seg(:,1),ref_seg(2,:), 'o-')
% plot(t_seg(:,1),ref_seg(3,:), 'o-')
% plot(t_seg(:,1),ref_seg(4,:), 'o-')
% xlabel('tempo (s)'), ylabel('amp (pu)')
% subplot(3,1,3), hold on, grid on
% plot(t_seg(:,1),ref_seg(5,:), 'o-')
% plot(t_seg(:,1),ref_seg(6,:), 'o-')
% plot(t_seg(:,1),ref_seg(7,:), 'o-')
% xlabel('tempo (s)'), ylabel('fase (rad)')

%--------------------------------------------------------------------------------------------------------------------
% ESTIMADOR FASORIAL DFT
[Sa, Sb, Sc, Aa, Ab, Ac, phia, phib, phic] = estimador_dft(Va_seg, Vb_seg, Vc_seg, Fs, f0);

granN = (Fs/N);
gradeN = 0:granN:granN*(N-1);

% % Plot componente fundamental
% figure
% subplot(2,1,1), hold on, grid on
% plot(Aa(:,2), 'o-')
% plot(Ab(:,2), 'o-')
% plot(Ac(:,2), 'o-')
% title('Amplitude (pu) - DFT')
% ylabel('pu'), xlabel('Ciclo')
% 
% subplot(2,1,2), hold on, grid on
% plot(phia(:,2), 'o-')
% plot(phib(:,2), 'o-')
% plot(phic(:,2), 'o-')
% title('Fase (rad) - DFT')
% ylabel('rad'), xlabel('Ciclo')

%--------------------------------------------------------------------------------------------------------------------
% CÁLCULO DE COMPONENTES SIMÉTRICAS

bin = round(f0/granN)+1;       % local da componente fundamental no espectro

% COM DADOS FFT
[A0, A1, A2, phi0, phi1, phi2] = comp_simetricas(Aa(:,bin)', Ab(:,bin)', Ac(:,bin)', phia(:,bin)', phib(:,bin)', phic(:,bin)');

% COM DADOS DE REFERÊNCIA (para teste do estimador)
[A0r, A1r, A2r, phi0r, phi1r, phi2r] = comp_simetricas(ref_seg(2,:), ref_seg(3,:), ref_seg(4,:),ref_seg(5,:), ref_seg(6,:), ref_seg(7,:));


% % Plot fasores da componente simetricas
% figure
% subplot(2,3,1), hold on, grid on
% plot(A0r, 's-'), plot(A0, 'o-')
% title('Amplitudes (pu)'), ylabel('pu'), xlabel('Ciclo')
% subplot(2,3,4), hold on, grid on
% plot(rad2deg(phi0r), 's-'),plot(rad2deg(phi0), 'o-')
% title('Fase (º)'), ylabel('º'), xlabel('Ciclo')
% subplot(2,3,2), hold on, grid on
% plot(A1r, 's-'), plot(A1, 'o-')
% title('Amplitudes (pu)'), ylabel('pu'), xlabel('Ciclo')
% subplot(2,3,5), hold on, grid on
% plot(rad2deg(phi1r), 's-'), plot(rad2deg(phi1), 'o-')
% title('Fase (º)'), ylabel('º'), xlabel('Ciclo')
% subplot(2,3,3), hold on, grid on
% plot(A2r, 's-'), plot(A2, 'o-')
% title('Amplitudes (pu)'), ylabel('pu'), xlabel('Ciclo')
% subplot(2,3,6), hold on, grid on
% plot(rad2deg(phi2r), 's-'), plot(rad2deg(phi2), 'o-')
% title('Fase (º)'), ylabel('º'), xlabel('Ciclo')

%--------------------------------------------------------------------------------------------------------------------
% ESTIMADOR DE FREQUÊNCIA - Derivada do ângulo do fasor de seq. positiva

% COM DADOS FFT
[freq_final_fft] = estimador_freq_deriv_ang(phi0, phi1, phi2, f0, Fs, ref_seg, m, N);

% figure, hold on, grid on
% plot(ref_seg(1,1:end), 'o-')
% plot(freq_final_fft, 'o-')
% title('Frequência (Hz)'), legend('Referência', 'Estimativa (DFT)');
% ylabel('Hz'), xlabel('Ciclo')



%--------------------------------------------------------------------------------------------------------------------
% CÁLCULO DE ERROS

% Seleção automática da componente fundamental
[num_wins,win_size] = size(Aa);
% gran = Fs/win_size/f0;       % granularidade
% bin = round(1/gran)+1;       % local da componente fundamental no espectro

% define indices (recorte dos dados)
i=1;  % inicio
f=length(ref_seg(1,:));%floor(T/(1/f0))-1; % fim (numero de ciclos)

% Armazena estimações em novos vetores (para vários estimadores)
% amplitudes
amp_med(1,i:f) = Aa(i:f,bin);  % dft

% fases
fase_med(1,i:f) = phia(i:f,bin);   

% frequência
freq_med(1,i:f) = freq_final_fft(1,i:f);


[tve, amp_error, phase_error, freq_error] = calcula_erros(amp_med, fase_med, freq_med, ref_seg(:,i:f));




% % % Plot erro frequência
% figure
% subplot(2,1,1),hold on, grid on  
% plot(ref_seg(1,:), 'o-')
% plot(freq_med(1,i:f), 'o-')
% title('Frequência (Hz)'), ylabel('Hz'), xlabel('Ciclo')
% subplot(2,1,2),hold on, grid on  
% plot(freq_error(1,:), 'o-')
% title('Erro de Frequência (Hz)'),% ylabel('Hz'), xlabel('Ciclo')
% 
% % Plot fasores da componente fundamental
% figure
% subplot(2,1,1), hold on, grid on
% plot(t_seg(i:f,1),ref_seg(2,:))%,'o-')
% plot(t_seg(i:f,1), amp_med(1,:))%, '*-')
% ylabel('amplitude (pu)'), xlabel('time (s)')
% legend('ref', 'DFT')
% title(strcat(['Errors, with f = ',num2str(freqs(1)),' Hz']))
% subplot(2,1,2), hold on, grid on
% plot(t_seg(i:f,1),rad2deg(ref_seg(5,:)))%, 'o-')
% plot(t_seg(i:f,1),rad2deg(fase_med(1,:)))%, '*-')
% ylabel('phase (rad)'), xlabel('time (s)')
% legend('ref', 'DFT')
% 
% 
% % Plot ERRO de fasores 
% figure  
% subplot(3,1,1), hold on, grid on%
% plot(t_seg(i:f,1), tve(1,:), '*-')
% ylabel('TVE (%)'), xlabel('tempo (s)')
% subplot(3,1,2), hold on, grid on
% plot(t_seg(i:f,1),amp_error(1,:), '*-')
% ylabel('Erro de amplitude (%)'), xlabel('tempo (s)')
% subplot(3,1,3), hold on, grid on
% plot(t_seg(i:f,1),phase_error(1,:), '*-')
% ylabel('Erro de fase (º)'), xlabel('tempo (s)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1
   max_tve(i) = max(tve(i,:));
   max_amp_error(i) = max(amp_error(i,:)); 
   max_phi_error(i) = max(phase_error(i,:)); 
end
%
nfasors = length(amp_error)
ymag = 0.025*vmrand(0,2,[1,nfasors])/pi;

yang = 0.025*vmrand(0,5,[1,nfasors])/pi;

yfreq = 0.0015*vmrand(0,2,[1,nfasors])/pi;

% Normaliza freq de casos

% % Aplica MLE etimator
% mle(y)
% mle(amp_error(1,:))

figure, hold on, grid on
histogram(yfreq, 'NumBins',100,'BinWidth',0.00005, 'Normalization', 'probability')
histogram(freq_error(1,:),'NumBins',100,'BinWidth',0.00005, 'Normalization', 'probability')

figure
subplot(1,2,1), hold on, grid on
histogram(ymag, 'NumBins',100,'BinWidth',0.001, 'Normalization', 'probability')
histogram(amp_error(1,:),'NumBins',100,'BinWidth',0.001, 'Normalization', 'probability')
title(strcat(['magnitude dispersion, with f = ',num2str(freqs(1)),' Hz'])), ylabel('probability'), xlabel('error (%)')
legend('VM\kappa=2', 'DFT')
subplot(1,2,2), hold on, grid on
histogram(yang, 'NumBins',100,'BinWidth',0.001, 'Normalization', 'probability')
histogram(phase_error(1,:),'NumBins',100,'BinWidth',0.001, 'Normalization', 'probability')
title(strcat(['phase dispersion, with f = ',num2str(freqs(1)),' Hz'])), ylabel('probability'), xlabel('error (rad)')
legend('VM\kappa=5', 'DFT')
