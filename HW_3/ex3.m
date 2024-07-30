
% submitted the mask that was drawn with the help of fiji photo
mask_path = "C:\Users\user\OneDrive - Bar-Ilan University - Students\Documents\Neurophotonics\Mask.tiff";
mask_tiff = Tiff(mask_path, 'r');
mask = read(mask_tiff);

%% q3
% find the read noise 
read_noise_path = 'D:\talia_chemda\read_noise_Gain24_expT0.021ms';
read_noise = Tiff2Matrix(read_noise_path);
read_noise_matrix = stdfilt(read_noise, true(7, 7)).^2;
mean_read_out_noise = mean(read_noise_matrix,3);

%% q4
% calculate the pixels non uniformity - we used the laser regular recording
% and applied it for both breath and regular recordings in future analysis
background_path = 'D:\talia_chemda\background_gain24_expt15ms_chemda';
path = 'D:\talia_chemda\fr20hz_gain24_expt15ms_chemda';
background_rec = Tiff2Matrix_2(background_path);
rec = Tiff2Matrix_2(path);
meanFrameBackground = mean(background_rec, 3);
meanFrame = mean(rec, 3) - meanFrameBackground;
windowSize = 7;
pixels_non_uniformity = stdfilt(meanFrame, true(windowSize)).^2; 


%% q5
% Inputs
N = 12;
saturation_capacity = 10500; 
g_in_dB = 24; 
% calculate G
G_base = 2^N / saturation_capacity;
G_in = 10^(g_in_dB / 20);
G = G_base * G_in;
disp(['G [DU/e]: ', num2str(G)]);


%% q6 + q7 + q8
Kq = 1 / (12);
% the short recording analysis
files_path_laser = 'D:\talia_chemda\fr20hz_gain24_expt15ms_chemda';
files_laser = dir(strcat(files_path_laser, '\*.tiff'));
for i = 1:length(files_laser)
    t = Tiff(strcat(files_path_laser, '\', files_laser(i).name), 'r');
    imageDataLaser = read(t);
    % subtract background from each frame
    imageDataLaser = imageDataLaser - uint16(meanFrameBackground);
    mean_window = double(imboxfilt(imageDataLaser, 7));
    mean_window = mean(nonzeros(mean_window));
    % calculate the raw contrast per window
    raw = (stdfilt(imageDataLaser, true(7, 7)));
    % calculate fixed contrast per window
    Kf = (raw - (mean_window .* G) - mean_read_out_noise - pixels_non_uniformity - Kq) ./ (mean_window.^2); 
    % extract only the data within the ROI - mask area
    Kf = Kf .* double(repmat(mask,[1 1 size(imageDataLaser, 3)]));
    Kf_matrix_laser(i) = mean(nonzeros(Kf(:)));
end

% the breath recording analysis
files_path_breath = 'D:\talia_chemda\fr20hz_gain24_expt15ms_chemda_breath';
files_breath = dir(strcat(files_path_breath, '\*.tiff'));
for i = 1:length(files_breath)
    t = Tiff(strcat(files_path_breath, '\', files(i).name), 'r');
    imageData = read(t);
    % subtract background from each frame
    imageData = imageData - uint16(meanFrameBackground);
    mean_window = double(imboxfilt(imageData, 7));
    mean_window = mean(nonzeros(mean_window));
    % calculate the raw contrast per window
    raw = (stdfilt(imageData, true(7, 7)));
    % calculate fixed contrast per window
    Kf = (raw - (mean_window .* G) - mean_read_out_noise - pixels_non_uniformity - Kq) ./ (mean_window.^2); 
    % extract only the data within the ROI - mask area
    Kf = Kf .* double(repmat(mask,[1 1 size(imageData, 3)]));
    Kf_matrix_breath(i) = mean(nonzeros(Kf(:)));
end

% plot for the short recording
time_per_frame = 1/20;
% convert to seconds for x axis
time_seconds = [1:length(Kf_matrix_laser)] * time_per_frame;
plot(time_seconds, Kf_matrix_laser, 'b-', 'LineWidth', 1.5); 
xlabel('Seconds'); 
ylabel('K_f^2'); 
title('Fixed Contrast K_f^2 Over Time - Short Recording'); 

% plot for breath recording
time_per_frame = 1/20;
% convert to seconds for x axis
time_seconds = [1:length(Kf_matrix_breath)] * time_per_frame;
plot(time_seconds, Kf_matrix_breath, 'b-', 'LineWidth', 1.5); 
% label when it appears in the plot the breath hold started and ended
% (specifically did not use the times we calculated because I wanted to
% show it from the plot)
xline([76.3,135.35, 197.25, 252.1],'--r','Breath hold start');
xline([119.1,176.25, 234.6],'--g','Breath hold end');
xlabel('Seconds'); 
ylabel('K_f^2'); 
title('Fixed Contrast K_f^2 Over Time - Breath Recording'); 
