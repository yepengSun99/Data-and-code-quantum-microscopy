clear; clc; close all;

%% ================= 1. Parameter Settings ================= %%
% TDMS file
file_quantum_classical = "classical.tdms";

Image_width  = 500;      % Image width: 500 pixels
Image_height = 400;      % Image height: 400 pixels

Dot_size = Image_width;  % Number of pixels per row
Capture_Mag = 2000;      % Samples acquired per pixel

fs = 2e6;               % Sampling frequency: 2 MHz
f_target = 500e3;        % Carrier frequency: 500 kHz

remove_dc = true;       % Enable pixel-wise DC offset removal

num_lines = Image_height;
N = Capture_Mag;

% Number of samples per image row
line_len = Dot_size * N;

%% ================= 2. FFT Frequency Bin Selection ================= %%
f_axis = (0:N-1) * fs / N;

[~, f_idx] = min(abs(f_axis - f_target));
f_actual = f_axis(f_idx);

fprintf('\n================ FFT Parameters ================\n');
fprintf('Image size: %d x %d pixels (Width x Height)\n', ...
    Image_width, Image_height);
fprintf('Image matrix size: %d x %d\n', ...
    Image_height, Image_width);
fprintf('Samples per pixel: %d\n', N);
fprintf('Sampling frequency: %.3f MHz\n', fs/1e6);
fprintf('Frequency resolution: %.3f kHz\n', fs/N/1e3);
fprintf('Target frequency: %.3f kHz\n', f_target/1e3);
fprintf('Selected frequency: %.3f kHz\n', f_actual/1e3);
fprintf('MATLAB FFT array index: %d\n', f_idx);
fprintf('DC removal enabled: %d\n', remove_dc);
fprintf('Scan order: odd rows forward, even rows reverse\n');

%% ================= 3. Read TDMS Data ================= %%
fprintf('\nReading TDMS data...\n');

signal = read_tdms_single_channel(file_quantum_classical);

expected_samples = num_lines * line_len;

fprintf('Actual number of samples: %d\n', numel(signal));
fprintf('Expected number of samples: %d\n', expected_samples);

if numel(signal) < expected_samples
    error(['Insufficient TDMS data: %d samples available, ' ...
           'but %d samples required.'], ...
           numel(signal), expected_samples);
end

if numel(signal) > expected_samples
    fprintf('Additional samples detected. Using the first %d samples.\n', ...
        expected_samples);
end

%% ================= 4. Pixel-wise FFT Reconstruction ================= %%
fprintf('\nStarting 500 kHz FFT image reconstruction...\n');

Img_500kHz_raw = reconstruct_500kHz( ...
    signal, num_lines, line_len, Dot_size, N, f_idx, remove_dc);

fprintf('500 kHz FFT image reconstruction completed.\n');

%% ================= 5. Colormap Configuration ================= %%
if exist('slanCM', 'file') == 2
    cmap = flipud(slanCM('amp', 256));
else
    cmap = gray(256);
end

%% ================= 6. Image Visualization ================= %%
figure('Color', 'w', ...
    'Name', '500 kHz Bidirectional FFT Reconstruction', ...
    'Position', [100 100 800 650]);

imagesc(Img_500kHz_raw);

axis image;
colormap(cmap);
colorbar;

xlabel('Column (pixel)');
ylabel('Row (pixel)');

title('Reconstructed 500 kHz FFT Magnitude');

%% ================= 7. Save Reconstructed Image ================= %%
save_name = 'Classical_500kHz_raw.mat';

save(save_name, 'Img_500kHz_raw');

fprintf('\n================ Processing Completed ================\n');
fprintf('Image matrix size: %d x %d\n', ...
    size(Img_500kHz_raw, 1), size(Img_500kHz_raw, 2));
fprintf('Even rows have been reversed for bidirectional scanning.\n');
fprintf('Reconstructed image saved to: %s\n', save_name);

plotedit on;

clear signal;

%% ===================================================================== %%
%% Function 1: Read a Single Channel from a TDMS File
%% ===================================================================== %%
function signal = read_tdms_single_channel(filepath)

    data = tdmsread(filepath);

    if iscell(data)
        table_data = data{1};
    elseif istable(data)
        table_data = data;
    else
        error('Unsupported TDMS data format. Expected a table or cell array.');
    end

    variable_name = table_data.Properties.VariableNames{1};

    signal = double(table_data.(variable_name));
    signal = signal(:);

    fprintf('Selected channel: %s\n', variable_name);
    fprintf('Total samples loaded: %d\n', numel(signal));
    fprintf('Signal data type: %s\n', class(signal));

end

%% ===================================================================== %%
%% Function 2: Pixel-wise FFT Reconstruction at 500 kHz
%% Odd rows: forward acquisition, retain pixel order
%% Even rows: reverse acquisition, reverse pixel order
%% ===================================================================== %%
function Img = reconstruct_500kHz( ...
    signal, num_lines, line_len, Dot_size, N, f_idx, remove_dc)

    Img = zeros(num_lines, Dot_size);

    signal = signal(:);
    total_samples = numel(signal);

    for line_index = 1:num_lines

        %% Step 1: Extract temporal samples for the current row
        line_start = (line_index - 1) * line_len + 1;
        line_end = line_start + line_len - 1;

        if line_end > total_samples
            error('Insufficient data for row %d.', line_index);
        end

        line_data = signal(line_start:line_end);

        %% Step 2: Reshape into individual pixel signals
        % Each column contains N temporal samples for one pixel.
        X = reshape(line_data, N, Dot_size);

        %% Step 3: Remove DC offset independently for each pixel
        if remove_dc
            X = X - mean(X, 1);
        end

        %% Step 4: Compute FFT magnitude spectrum
        Y = abs(fft(X, [], 1));

        %% Step 5: Extract the 500 kHz FFT magnitude
        row_values = Y(f_idx, :);

        %% Step 6: Correct bidirectional scanning order
        % Rows 1, 3, 5, ...: retain original order.
        % Rows 2, 4, 6, ...: reverse pixel order.
        if mod(line_index, 2) == 0
            row_values = fliplr(row_values);
        end

        Img(line_index, :) = row_values;

        %% Step 7: Display reconstruction progress
        if mod(line_index, 50) == 0 || line_index == num_lines
            fprintf('Reconstruction progress: %d / %d rows\n', ...
                line_index, num_lines);
        end

    end

end