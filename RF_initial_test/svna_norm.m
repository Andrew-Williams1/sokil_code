% despriction of what it does

% output = [freq_avg;mag_avg;phase_avg;re;im]
clear all

% x = indexing "CSV_.csv" file from 10-21 test
x = [16 17 13 15];

% 11 = empty chamber w/ tube
normData = svna_data_analysis(11);

% Set up Legend
legend_all = ["Empty chamber (w/ cardboard)", "Empty chamber (w/out cardboard)", "Steel ball", "Acrylic ball", "HDPE Green sheet", "Steel block", "Aluminum sheet"];
% all possible x values
x_all = [11 12 13 15 17 16 18];
legend_L = []; % Legend string list
for y = x
    i = find(x_all==y); % index of desired x
    legend_L = [legend_L, legend_all(i)];
end
%   
% Make plots
linewidth = .9;
figure
for i = x
    data = svna_data_analysis(i);
    freq = data(1,:);
    % comp = complex
    comp = (data(4,:)+j*data(5,:))./(normData(4,:)+j*normData(5,:));
    re = real(comp);
    im = imag(comp);
    mag = sqrt(re.^2 + im.^2);
    logmag = 20*log10(mag); % log magnitude
    phase = atan2(im,re);
    
    subplot(3,1,1);
    hold on;
    plot(freq, logmag, 'LineWidth', linewidth);
    xlabel('Frequency [MHz]');
    ylabel('Normalized Log Magnitude');
    title('Normalized Log Magnitude |S11| comparison');
    legend(legend_L, 'Location', 'southwest');

    hold off;
    subplot(3,1,2);
    hold on;
    plot(freq, unwrap(phase), 'LineWidth', linewidth);
    xlabel('Frequency [MHz]');
    ylabel('Normalized Unwrapped Phase [rads]');
    title('Normalized S11 Phase comparison');
    legend(legend_L, 'Location', 'southwest');
    hold off;

    subplot(3,1,3);
    hold on;
    fdiff = (10^7)/3; % svna from 700 MHz - 3 GHz w/ 750 data points
    N = 2000; % N > 2*750
    fs = N*fdiff;
    dt = 1/fs;
    t = 0:dt:(N)*dt;
    f = -fs/2:fs/N:fs/2;
    mask = (abs(f) > 700*10^6) .* (abs(f) < 3*10^9); % redundant
    FD = mask.*zeros(1,length(f));
    %comp = [zeros(1,150), comp];
    vna_start = 1152; % index of 503.33 MHz
    vna_end = 1901; % index of 3 GHz
    vna_negative_start = 101; % -3
    vna_negative_end = 850; % -503.3
    FD(vna_start:vna_end) = comp;  % insert the 750 real and imag samples into the positive frequencies
    FD(vna_negative_start:vna_negative_end) = flip(conj(comp));  % do the same to negative freqs.
    TD = ifft(FD); % normalized time domain signals
    plot(t,real(TD), 'LineWidth', linewidth);
    %hold on
    %plot(t,imag(TD), 'LineWidth', linewidth);
    xlabel('Wrong units');
    ylabel('Normalized S11');
    title('Normalized Inverse Fourier Transform');
    legend(legend_L, 'Location', 'northeast');
    hold off;  
end

    

    

