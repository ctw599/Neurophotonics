% This code takes as input the recordings given with tiff files.
% It returns the temporal noise, global spatial noise, local spatial noise
% averged in time and local spatial noise averaged over all frames.
% Input: recording paths
% Output: temporalNoise, globalSpatialNoise, localSpatialNoise_c,
% localSpatialNoise_d, temporalNoise2_Mean, theoreticalValue

% Paths to the recording folders
recordingPaths = {'C:\Users\user\Downloads\records\records\1. WithCover_Gain0dB_expT0.021ms_BlackLevel30DU', ...
                  'C:\Users\user\Downloads\records\records\2. WithCover_Gain24dB_expT0.021ms_BlackLevel30DU', ...
                  'C:\Users\user\Downloads\records\records\3. WithCover_Gain24dB_expT10ms_BlackLevel30DU', ...
                  'C:\Users\user\Downloads\records\records\4. WhitePaper_Gain24dB_expT0.5ms_BlackLevel0DU'};
temporalNoise = [];
globalSpatialNoise = [];
localSpatialNoise_c = [];
localSpatialNoise_d = [];
validRecordingPaths = {};
% Loop through each recording path
for i = 1:length(recordingPaths)
    folderPath = recordingPaths{i};
    % Load frames using provided function
    rec = Tiff2Matrix(folderPath);  
    
    if isempty(rec)
        warning('No frames loaded from folder: %s', folderPath);
        continue;
    else
        % Add to valid paths - just to check variables
        validRecordingPaths = [validRecordingPaths; folderPath]; 
    end
    
    numFrames = size(rec, 3);
    [rows, cols] = size(rec(:,:,1));
    
    % Temporal Noise
    % Standard deviation across the third dimension (time) - as written in
    % the powerpoint
    temporalNoise_rec = std(rec, 0, 3); 
    temporalNoise = [temporalNoise; mean(temporalNoise_rec(:))];
    
    % Global Spatial Noise (after averaging in time)
    meanFrame = mean(rec, 3);
    globalSpatialNoise = [globalSpatialNoise; std(meanFrame(:))];
    
    % Local Spatial Noise with window size of 7 (after averaging in time)
    windowSize = 7;
    localSpatialNoise_c_rec = stdfilt(meanFrame, true(windowSize)); 
    localSpatialNoise_c = [localSpatialNoise_c; mean(localSpatialNoise_c_rec(:))];
    
    % Local Spatial Noise with window size of 7 per frame - then average over all frames
    localSpatialNoise_d_rec = zeros(rows, cols, numFrames);
    
    for j = 1:numFrames
        localSpatialNoise_d_rec(:,:,j) = stdfilt(rec(:,:,j), true(windowSize)); 
    end
    localSpatialNoise_d = [localSpatialNoise_d; mean(localSpatialNoise_d_rec(:))];
end
% Calculate total noise
totalNoise = sqrt(temporalNoise.^2 + localSpatialNoise_d.^2);
% Display results
tableData = table(validRecordingPaths, temporalNoise, globalSpatialNoise, localSpatialNoise_c, localSpatialNoise_d, totalNoise, ...
                  'VariableNames', {'RecordingPath', 'TemporalNoise', 'GlobalSpatialNoise', 'LocalSpatialNoise_c', 'LocalSpatialNoise_d', 'TotalNoise'});
disp(tableData);


% Path to the specific recording folder - the last recording with the
% camera opened
folderPath = recordingPaths{4};

% Load frames using provided function
rec = Tiff2Matrix(folderPath);

if isempty(rec)
    error('No frames loaded from folder: %s', folderPath);
end

numFrames = size(rec, 3);

% Temporal Noise
% Standard deviation across the third dimension (time)
temporalNoise2 = temporalNoise(4);
meanSignal = mean(rec, 3);

% Calculate the ratio of temporal noise squared to mean signal
temporalNoise2_Mean = (temporalNoise2 .^ 2) ./ meanSignal;
% Average over all pixels
temporalNoise2_Mean = mean(temporalNoise2_Mean(:)); 

% Given parameters - as stated in the homework file
nBits = 12; 
maxCapacity = 10500; 
gain_dB = 24; 

% Theoretical value calculation
theoreticalValue = (2^nBits / maxCapacity) * 10^(gain_dB / 20);

fprintf('Calculated temporalNoise / Mean: %.5f\n', temporalNoise2_Mean);
fprintf('Theoretical value: %.5f\n', theoreticalValue);

% Calculate the mean image for the fourth recording
meanImage = mean(rec, 3);
figure;
imagesc(meanImage);
axis image;
colormap('gray');
title('Mean Image');


