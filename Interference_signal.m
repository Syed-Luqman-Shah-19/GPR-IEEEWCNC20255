%% This repository contains simulation code for the paper titled 'Interference Prediction Using Gaussian Process Regression and Management Framework for Critical Services in % Local 6G Networks'. The paper was authored by Syed Luqman Shah, Nurul Huda Mahmood, and Matti Latva-aho. The paper has been accepted for publication in IEEE WCNC 2025 % Conference.
 % email: sayedluqmans@gmail.com
%%%%%%%%%%%%%%
% This code generates correlated interference signal. Where each value of
% the current sample impacted by the past 'filterLength' values.
% The correlation profile is given by 'filterCoefs'
% Here a simple uniform filter is applied. It can be changed to other
% filter types, e.g. the Bessel function, Jakes model, etc.


%%
% function [sum_interference,imf,residual] = correlated_signal_func()
clc;
clearvars;
close all;

%nrOfSamples = 1e4; % number of samples to generate
nrOfSamples = 200; % number of samples to generate
inr_db = [5, 2, 0, -3, -10, 1]; % mean power (normalized by noise) of the each 
                             % interferer in dB
inr_lin = 10.^(inr_db/10);   % cnvert to linear
nrOfInterferers = length(inr_db);  % number of interferers

filterLength = 40;
filterCoefs = 1/filterLength*ones(filterLength,1);
filterStateI = zeros(filterLength - 1, nrOfInterferers);
filterStateQ = zeros(filterLength - 1, nrOfInterferers);

rng default

samplingFreq = 4e6;
sampleTime = 1/samplingFreq;

time = sampleTime*(0:nrOfSamples-1);
sum_interference = 0; % the sum interference signal

for ii = 1:nrOfInterferers
    noiseI = randn(nrOfSamples,1); % normally distributed inPhase component
    noiseQ = randn(nrOfSamples,1); % normally distributed Quadrature component
    
    [noiseI, stateI] = filter(filterCoefs, 1, noiseI, filterStateI(:,ii));
    [noiseQ, stateQ] = filter(filterCoefs, 1, noiseQ, filterStateQ(:,ii));
    filterStateI(:, ii) = stateI;
    filterStateQ(:, ii) = stateQ;
    
    nn = noiseI + 1i*noiseQ; % the overall correlated signal
    interfSignal = (real(nn).^2 + imag(nn).^2); % the power of the signal
    interfSignal = interfSignal/mean(interfSignal)*inr_lin(ii);
        % normalize the power and multiply by the mean interference power
        % (normalized by the noise power)
    
    sum_interference = sum_interference + interfSignal;
end

[imf,residual] = emd(sum_interference);
csvwrite('interference.csv', sum_interference);
% end


%%% FOR COMPARISON - An uncorrelated signal

% sum_interfUncorr = 0;
% 
% for ii = 1:nrOfInterferers
%     noiseI = randn(nrOfSamples,1);
%     noiseQ = randn(nrOfSamples,1);
% 
%     nn = noiseI + 1i*noiseQ;
%     interfSignal = (real(nn).^2 + imag(nn).^2);
%     interfSignal = interfSignal/mean(interfSignal)*inr_lin(ii);
%     
%     sum_interfUncorr = sum_interfUncorr + interfSignal;
% end


%% PLOT THE TWO SIGNALS
figure;
% plot(sum_interference(filterLength + 1:filterLength + min(length(sum_interference), 100)), 'b', 'LineWidth',2)
plot(sum_interference, 'b', 'LineWidth',2)
fontsize(gca, 16, 'points')
xlabel('Time (s)')
ylabel('Amplitude')
% hold on; grid on;
% plot(sum_interfUncorr(filterLength + 1:filterLength + min(length(sum_interfUncorr), 100)), 'r', 'LineWidth',2)
% legend('Correlated Signal', 'Uncorrelated Signal')


%% AUTOCRRELATION
[acf1, lags1] = autocorr(sum_interference, NumLags=199)
% [acf1, lags1] = autocorr(sum_interference)

figure
stem(lags1, acf1);
fontsize(gca, 16, 'points')
% title('Autocorrelation variation of the correlated interference signal')
ylabel('Autocorrelation', 'FontSize', 16)
xlabel('Time lag', 'FontSize', 16)

%% SIGNAL
clc;
clearvars;
nrOfPreSamples=200;
nrOfSamples = 200;

%inr_db = [5, 3, 0, -2, -5, 10, -10, 1, -1, 4, -3, -9]; % mean INR of the different interferers in dB >> you can change these values as you wish.can be in range of [-10,5]
inr_db = [20]; % generate signal with INR value of 20 dB
inr_lin = 10.^(inr_db/10); % converting from dB to linear
nrOfInterferers = length(inr_db);

interf_matrix = zeros(nrOfSamples, nrOfInterferers);
% generate the random interference signals with the given SNRs for each interferer. Here we assume independent Rayleigh fading. (The power of Rayleigh fading is exponentially distributed.)
for ii = 1:nrOfInterferers
%     interf_matrix(:,ii) = exprnd(inr_lin(ii),[nrOfPreSamples,1]);
    interf_matrix(:,ii) = randn([nrOfPreSamples,1]);
end
% sum the individual interference signals to get the sum interference signal.
signal = sum(interf_matrix, 2);
%[imf,residual] = emd(sum_interference)

figure;
% plot(sum_interference(filterLength + 1:filterLength + min(length(sum_interference), 100)), 'b', 'LineWidth',2)
plot(signal, 'b', 'LineWidth',2)
fontsize(gca, 16, 'points')
xlabel('Time (s)')
ylabel('Amplitude')
% hold on; grid on;



