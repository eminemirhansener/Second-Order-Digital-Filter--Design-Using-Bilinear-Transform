% Parameters
f_c = 500; % cutoff frequency
f_s = 2500; % sampling frequency
zeta = 0.7;

frequencies = [1, 450, 800, 1000, 1500]; % frequency of applied sinusoidal signal

A = 1;
phi = 0;
T = 1/f_s;
N = f_s * 1;

%% Calculating coefficients
w_c = 2 * pi * f_c;
w_c_prime = (2 / T) * tan(w_c * (T / 2));

a0 = (4 / (T^2)) + ((4 / T) * zeta * w_c_prime) + (w_c_prime^2);
a1 = 2 * (w_c_prime^2) - (8 / (T^2));
a1 = a1 / a0;
a2 = w_c_prime^2 + (4 / (T^2)) - ((4 / T) * zeta * w_c_prime);
a2 = a2 / a0;

b0 = (w_c_prime^2);
b0 = b0 / a0;
b1 = (2 * (w_c_prime^2));
b1 = b1 / a0;
b2 = w_c_prime^2;
b2 = b2 / a0;

a0 = a0 / a0;

%% Input signal
n = 0:1:N-1;

figure

for i = 1:length(frequencies)
    f = frequencies(i);
    x = A * sin(2 * pi * f * n * T + phi);
    y = zeros(N, 1);
    
    for j = 3:N
        y(j) = b0 * x(j) + b1 * x(j-1) + b2 * x(j-2) - a1 * y(j-1) - a2 * y(j-2);
    end

    subplot(length(frequencies), 1, i);
    plot(n * T, x, 'b', n * T, y, 'r')
    title(sprintf('%d Hz input and Filtered Output Signals', f));
    xlabel('Time');
    ylabel('Amplitude');
    legend('Input Signal', 'Filtered Signal', 'Location', 'best');
    grid on;
end

a = [a0 a1 a2];
b = [b0 b1 b2];
[H, f] = freqz(b, a, [], f_s);

figure
y_cutoff = -3;
x_cutoff = interp1(20 * log10(abs(H)), f, y_cutoff, 'linear');

y_stopband = -40;
x_stopband = interp1(20 * log10(abs(H)), f, y_stopband, 'linear');

y_passband = -0.5;
x_passband = interp1(20 * log10(abs(H)), f, y_passband, 'linear');

disp('Cutoff Frequency (at -3 dB):');
disp(x_cutoff);

disp('Stopband Frequency (at -40 dB):');
disp(x_stopband);

disp('Passband Frequency (at -0.5 dB):');
disp(x_passband);

plot(f, 20 * log10(abs(H)));
title('Frequency-Magnitude Response of Discrete Filter Design');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;