
function [dHbR, dHbO, fig] = CalcNIRS(dataFile, SDS, tissueType, extinctionCoefficientsFile, DPFperTissueFile, relDPFfile, plotChannelIdx)
    % CalcNIRS - calculate and plot HbR HbO
    % Input:
    % dataFile - .mat file with intensity data.
    % DS.Lambda : two wavelengths (in nm)
    % t : time vector
    % d : intensity data of 20 channels
    % 20 first rows-> first wavelength, 20 last rows->second wavelength
    % SDS - Source-Detector Separation distance in cm
    % tissueType - one of the rows in DPFperTissueFile
    % plotChannelIdx - vector with numbers in the range of [1-20] indicating channels to plot
    % extinctionCoefficientsFile - .csv file with extinction coefficients data
    % DPFperTissueFile - .txt file with tissue-specific DPF values
    % relDPFfile - .csv file with relative DPF coefficients according to wavelength
    % Output :
    % dHbR - HbR concentration change for all channels (nx20) where n is time vector length
    % dHbO - HbO concentration change for all channels (nx20) where n is time vector length
    % fig - handle to figure. Empty if plotChannelIdx==[]

    if nargin < 7
        plotChannelIdx = [];
    end

    % Load the .mat data file
    data = load(dataFile);
    Lambda = data.SD.Lambda;
    timeVector = data.t;
    intensityData = data.d;
    
    % Read the extinction coefficients from the CSV file
    extinctionData = readtable(extinctionCoefficientsFile);
    extinctionData = renamevars(extinctionData,'x_SpecificExtinctionCoeffsTakenFromBORLWaterMeasuredByMatcherAt','Wavelength'); 
    extinctionData = renamevars(extinctionData,'Var2','Water'); 
    extinctionData = renamevars(extinctionData,'Var3','HbO2'); 
    extinctionData = renamevars(extinctionData,'Var4','HHb'); 
    extinctionData = renamevars(extinctionData,'Var5','FatSoybean'); 
    extCoeffLambda1 = extinctionData{extinctionData.Wavelength == Lambda(1), {'HbO2', 'HHb'}};
    extCoeffLambda2 = extinctionData{extinctionData.Wavelength == Lambda(2), {'HbO2', 'HHb'}};

    % Read the DPF data for the given tissue type
    DPFData = readtable(DPFperTissueFile, 'Delimiter', '\t');
    DPF = DPFData{strcmp(DPFData.Tissue, tissueType), 'DPF'};

    % Read the relative DPF coefficients
    relDPFData = readtable(relDPFfile);
    relDPF1 = relDPFData{relDPFData.wavelength == Lambda(1), 'relDPFcoeff'};
    relDPF2 = relDPFData{relDPFData.wavelength == Lambda(2), 'relDPFcoeff'};

    % Calculate the actual DPF for each wavelength
    actualDPF1 = DPF * relDPF1;
    actualDPF2 = DPF * relDPF2;

    % Split intensity data into two wavelengths
    I1 = intensityData(:, 1:20); % First 20 columns -> first wavelength
    I2 = intensityData(:, 21:40); % Last 20 columns -> second wavelength

    % Convert intensities to optical densities
    OD1 = log(I1);
    OD2 = log(I2);

    % Calculate delta OD by subtracting the baseline (first measurement)
    dOD1 = [-OD1 + OD1(1,:)];
    dOD2 = [-OD2 + OD2(1,:)];

    pathlength_1 = actualDPF1 * SDS;
    pathlength_2 = actualDPF2 * SDS;

    % Create the extinction coefficient matrix
    % ASK IF TO DIVIDE OR TO MULTIPLY PATHLENGTH AND WHY TWO DIFFERENT WHEN
    % THE SLIDE SHOWS ONLY ONE DIFFERENT
    E = [extCoeffLambda1(1)/pathlength_1, extCoeffLambda1(2)/pathlength_1; ...
         extCoeffLambda2(1)/pathlength_2, extCoeffLambda2(2)/pathlength_2];

    % Inverse of the extinction coefficient matrix
    invE = inv(E);

    % Calculate changes in HbO and HbR concentrations for all channels
    dHbO = zeros(size(dOD1));
    dHbR = zeros(size(dOD1));

    for i = 1:20
        dOD = [dOD1(:,i), dOD2(:,i)];
        dC = invE * dOD';
        dHbO(:, i) = dC(1,:);
        dHbR(:, i) = dC(2,:);
    end

    % Plot if required
    if ~isempty(plotChannelIdx)
        fig = figure;
        hold on;
        for i = plotChannelIdx
            subplot(length(plotChannelIdx), 1, find(plotChannelIdx == i));
            plot(timeVector, dHbO(:,i), 'r', 'DisplayName', 'HbO');
            hold on;
            plot(timeVector, dHbR(:,i), 'b', 'DisplayName', 'HbR');
            title(['Channel ' num2str(i)]);
            xlabel('Time (s)');
            ylabel('Concentration Change');
            legend show;
        end
        hold off;
    else
        fig = [];
    end
    figure
    Ts = 2823/361.2160;
    y = fft(intensityData(:,1));
    fs = 1/Ts;
    f = (0:length(y)-1)*fs/length(y);
    plot(f,abs(y))
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title('Magnitude')
end

