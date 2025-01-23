%% This repository contains simulation code for the paper titled 'Interference Prediction Using Gaussian Process Regression and Management Framework for Critical Services in % % Local 6G Networks'. The paper was authored by Syed Luqman Shah, Nurul Huda Mahmood, and Matti Latva-aho. The paper has been accepted for publication in IEEE WCNC 2025 % % %%%%% Conference.
 % email: sayedluqmans@gmail.com

function [signal] = signal_only()
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
        interf_matrix(:,ii) = exprnd(inr_lin(ii),[nrOfPreSamples,1]);
    end
    % sum the individual interference signals to get the sum interference signal.
    signal = sum(interf_matrix, 2);
    %[imf,residual] = emd(sum_interference)
    csvwrite('signal.csv', signal);

end