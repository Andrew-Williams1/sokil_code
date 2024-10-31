% Sokil Clinic 2024-2025
% odshgjrpowejfgwoerisdlkjvs
% chat gpt
function create_plot_with_error(freq_avg, mag_avg, phase_avg, mag_SER, phase_SER);
    % Plot Magnitude
    figure;
    subplot(2,1,1);
    plot(freq_avg, mag_avg, 'DisplayName', 'Data');
    hold on;
    % Add error bars showing the standard error of the regression
    %errorbar(freq_avg, mag_avg, mag_SER .* ones(size(freq_avg)), 'r.', 'DisplayName', 'SER Error Bars');
    xlabel('Frequency [MHz]');
    ylabel('Magnitude |S11|');
    legend('Location', 'best');
    title('Magnitude |S11| vs. Frequency [MHz] with Standard Error of Regression');
    hold off;
    % Plot Phase
    subplot(2,1,2);
    plot(freq_avg, phase_avg, 'DisplayName', 'Data');
    hold on;
    % Add error bars showing the standard error of the regression
    errorbar(freq_avg, phase_avg, phase_SER .* ones(size(freq_avg)), 'r.', 'DisplayName', 'SER Error Bars');
    xlabel('Frequency [MHz]');
    ylabel('Phase [rads]');
    legend('Location', 'best');
    title('Phase vs. Frequency [MHz] with Standard Error of Regression');
    hold off;

    % Plot complex vectors
    complex = [10.^(mag_avg/20) .* exp(1j*phase_avg)];
    figure;
    %plot(x, y, 'DisplayName', 'Data');
    % Set plot labels and legend
    plot3(real(complex),imag(complex), freq_avg);
    xlabel('Real');
    ylabel('Imaginary');
    zlabel('Frequency [MHz]')
    title('Real vs. Imaginary');
    %scatter(1:807, freq_avg(2:808)-freq_avg(1:807));
end