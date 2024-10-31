% Sokil Clinic 2024-2025
% Takes .csv files to print data analysis plots
% RF plotter
% nanovna_data_analysis
% only works with 10-7 tests
% filename: sokil10-7_x.csv
% THREE COLUMNS:
% 1. Frequency (Hz)
% 2. Magnitude (linear or dB)
% 3. Phase (radians)

function output = nanovna_data_analysis(x)
    % x = list of indexes used in file naming convention
    % x = 1:5, 6:10, or 11:13
    for i = 1:length(x) % Goes through each whole file
        file_i = string(x(i));
        csv_filename = append('sokil10-7_',file_i,'.csv');
    
        % Read the CSV file into a matrix
        data = csvread(data_folder/csv_filename);
        
        % Extract the columns as matrices with all the data
        freq(i,:) = data(:, 1);  % frequency
        mag(i,:) = data(:, 2);  % magnitude 
        phase(i,:) = data(:, 3);  % phase
    end
    
    % Change frequency in MHz
    freq = freq .* (10^(-15));
        % Error regression: go through columns 2 and 3, combine align the 5 
        % values for each column, standard deviation for each point, find average
        % plot with error based on std
        
        % Fourier Transform: convet to complex numbers (keep as radians)
    
    % Calculate the standard deviations and averages at each point
    len = length(freq);
    for i = 1:len
        freq_std(i) = std(freq(:,i));
        freq_avg(i) = mean(freq(:,i));
        mag_std(i) = std(mag(:,i));
        mag_avg(i) = mean(mag(:,i));
        phase_std(i) = std(phase(:,i));
        phase_avg(i) = mean(phase(:,i));
    end
    freq_avg = freq_avg(1):((freq_avg(len)-freq_avg(1))/(len-1)):freq_avg(len);
    %mag_avg = abs(mag_avg);
    %mag_avg = db2mag(mag_avg);

    l = length(x); % set of data points at same frequency
    for n = 1:len % Goes thorugh each point
        % Calculate the residuals (difference between actual and mean values)
        residuals = mag(:,n)' - repmat(mag_avg(n),1,l);
        % Calculate the Standard Error of Regression (SER)
        mag_SER(n) = sqrt(sum(residuals.^2) / (len - 1));  % Standard error of regression
        % Calculate the residuals (difference between actual and mean values)
        residuals = phase(:,n)' - repmat(phase_avg(n),1,l);
        % Calculate the Standard Error of Regression (SER)
        phase_SER(n) = sqrt(sum(residuals.^2) / (len - 1));  % Standard error of regression
    end

    %create_plot_with_error(freq_avg, mag_avg, phase_avg, mag_SER, phase_SER);
    output = [freq_avg; mag_avg; phase_avg; mag_SER; phase_SER];
end
