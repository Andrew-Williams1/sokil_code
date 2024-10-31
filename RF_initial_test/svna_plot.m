
x = [13 15 17 16];

% Set up Legend
legend_all = ["Empty chamber (w/ cardboard)", "Empty chamber (w/out cardboard)", "Steel ball", "Acrylic ball", "HDPE Green sheet", "Steel block", "Aluminum sheet"];
% all possible x values
x_all = [11 12 13 15 17 16 18];
legend_L = []; % Legend string list
for y = x
    i = find(x_all==y); % index of desired x
    legend_L = [legend_L, legend_all(i)];
end

figure
for i = x
    data = svna_data_analysis(i);
    freq = data(1,:);
    mag = data(2,:);
    logmag = 20*log10(mag);
    phase = data(3,:); 
    % Log magnitude
    subplot(3,1,1);
    hold on;
    %plot(freq, logmag, 'LineWidth', 1.5);
    plot(freq, data(4,:), 'LineWidth', 1.5);
    xlabel('Frequency [MHz]');
    ylabel('Log Magnitude');
    title('Log Magnitude |S11| comparison');
    legend(legend_L, 'Location', 'southwest');
    hold off;

    % Phase
    subplot(3,1,2);
    hold on;
    %plot(freq, unwrap(phase), 'LineWidth', 1.5);
    plot(freq, data(5,:), 'LineWidth', 1.5);
    xlabel('Frequency [MHz]');
    ylabel('Unwrapped Phase [rap]');
    title('S11 Phase comparison');
    legend(legend_L, 'Location', 'southwest');
    hold off;

    % Inverse Fourier Transform
    subplot(3,1,3);
    hold on;
    % comp = complex
    comp = data(4,:) + j*data(5,:);
    for i = 1:length(freq)
        if freq(i) < 700
            comp(i) = 0;
        end
    end
    TD = ifft(comp); % time domain signal
    plot(abs(TD), 'LineWidth', 1.5);
    xlabel('Wrong units');
    ylabel('Normalized S11?');
    title('Inverse Fourier Transform Signal');
    legend(legend_L, 'Location', 'northeast');
    hold off;
end

