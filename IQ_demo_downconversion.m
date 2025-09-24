% Ideal_Downconversion_IQDemodulator.m
%
% This script simulates a downconversion and IQ demodulation
% process for a cavity signal using a custom-built FIR filter.
% The parameters used are similar to 650 MHz signal, the cavity signal is being phase modulated 
% with at 5 kHz. As of now the filters are creating issues for the IQ demodulated signal. 
clear  

% --- Parameters ---
tmax = .001; % Total simulation time in seconds.
Fs = 10e9; % Sampling frequency
dt = 1/Fs; % Time step
fcavity = 650e6; % Cavity resonant frequency
fLO1 = 645e6; % First local oscillator frequency
%IF_FREQ = fcavity - fLO1; % Intermediate frequency
t = 0:dt:tmax-dt; % Time vector
t=t';
w = 2*pi*[0:numel(t)/2-1 -numel(t)/2:-1].'/(numel(t)*mean(diff(t))) ;

% --- Signal Generation ---
% Corrected: Model the vibration as phase modulation.
cavity = sin(2*pi*fcavity*t + 0.2 * sin(2*pi*5000*t));

% Local Oscillator 1 (LO1)
LO1 = cos(2*pi*fLO1*t);

% --- Downconversion Stage 1: Mixer and Low-Pass Filter ---
% The mixer model is Mini-Circuits ZFM-2000+ which uses dioded
% The mixer diode model is shown in wikipedia article: https://en.wikipedia.org/wiki/Frequency_mixer
downconverted_signal =(cavity + LO1) + 0.5*(cavity.^2 + LO1.^2)+ cavity .* LO1;

%What does the FFT look like before the filter 
figure(1)
loglog(abs(w)/(2*pi*1e6),abs(fft(downconverted_signal)))
hold on 

% Design the first FIR low-pass filter
% An FIR filter order of 256 is a good starting point.
order = 256; 
first_LPF = 20e6;
cutoff_lp1_norm = first_LPF/ (Fs/2); 
b1 = my_fir_lowpass(order, cutoff_lp1_norm);

% Apply the first low-pass filter using the standard filter function
% Note: The FIR filter introduces a delay of (order/2) samples.
cavity_if_signal_raw = filter(b1, 1, downconverted_signal);

% Correct for the filter's delay
delay_samples = order / 2;
cavity_if_signal = [cavity_if_signal_raw(delay_samples+1:end,1); zeros(delay_samples,1)];
% for plot in figure(1)
loglog(abs(w)/(2*pi*1e6),abs(fft(cavity_if_signal)))
legend('Downconverted', [num2str(first_LPF/1e6) ' MHz LPF of Downconverted'])
legend('boxoff')


IF_FREQ=5e6;


% --- Downconversion Stage 2: IQ Demodulation ---------------------------
LO2_I = cos(2*pi*(IF_FREQ)*t);
LO2_Q = -sin(2*pi*(IF_FREQ)*t);
mixed_I = cavity_if_signal .* LO2_I;
mixed_Q = cavity_if_signal .* LO2_Q;

% Design the second low-pass filter for the baseband signals
order2=8192;
cutoff_demod_norm = (7000) / (Fs/2);
b2 = my_fir_lowpass(order2, cutoff_demod_norm);

% Apply the second low-pass filter
I_demod_raw = filter(b2, 1, mixed_I);
Q_demod_raw = filter(b2, 1, mixed_Q);

% Correct for the filter's delay
I_demod = [I_demod_raw(delay_samples+1:end); zeros(delay_samples,1)];
Q_demod = [Q_demod_raw(delay_samples+1:end); zeros(delay_samples,1)];

clear LO2_Q LO2_I cavity_if_signal_raw mixed_I mixed_Q 

% --- Plotting the Results ---
figure(2);
tiledlayout(3,1)

% Plot 1: Original RF signals
nexttile;
plot(t * 1e9, cavity);
title(sprintf('Original %.0f MHz Signal', fcavity/1e6));
xlabel('Time (ns)');
ylabel('Amplitude');
grid on;
xlim([0, 100]);

% Plot 2: Downconverted 5 MHz IF signals
nexttile;
plot(t * 1e6, cavity_if_signal);
title(sprintf('Downconverted %.0f MHz IF Signal', IF_FREQ/1e6));
xlabel('Time (µs)');
ylabel('Amplitude');
grid on;
xlim([0, 0.5]);

% Plot 3: Demodulated Baseband I and Q signals
nexttile;
plot(t * 1e3, I_demod, 'DisplayName', 'I-component');
hold on;
plot(t * 1e3, Q_demod, 'DisplayName', 'Q-component');
title('Demodulated Baseband I/Q Signals');
xlabel('Time (ms)');
ylabel('Amplitude');
grid on;
xlim([0, max(t)*1e3]); % Show 5 periods of the 1 kHz signal
legend('show');
hold off;


figure(3)
loglog(abs(w)/(2*pi*1e6),abs(fft(I_demod_raw)),'DisplayName', 'I-component')
hold on 
loglog(abs(w)/(2*pi*1e6),abs(fft(Q_demod_raw)),'DisplayName', 'Q-component')


% --- Custom FIR Low
% -Pass Filter Design (Without Toolbox) ---
% This function designs a low-pass FIR filter using the windowed-sinc method.
function b = my_fir_lowpass(order, cutoff_norm)
    % N: Filter order (number of coefficients)
    N = order;
    
    % Create the ideal sinc filter impulse response
    % Time vector for the filter's impulse response
    t_fir = -(N-1)/2 : (N-1)/2; 
    
    % Ideal impulse response (sinc function)
    sinc_impulse = sin(2 * pi * cutoff_norm * t_fir) ./ (pi * t_fir);
    
    % Handle the NaN at t=0
    sinc_impulse(t_fir == 0) = 2 * cutoff_norm;

    % Create a Hamming window to reduce ripple
    hamming_window = 0.54 - 0.46 * cos(2 * pi * (0:N-1) / (N-1));
    
    % Apply the window to the sinc function
    b = sinc_impulse .* hamming_window;

    % Normalize the filter coefficients so the DC gain is 1
    b = b / sum(b);
end