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

function [dHbR, dHbO, fig] = CalcNIRS(dataFile, SDS, tissueType, plotChannelIdx, extinctionCoefficientsFile, DPFperTissueFile, relDPFfile)
    if isfile(dataFile) && isfile(extinctionCoefficientsFile) && isfile(DPFperTissueFile) && isfile(relDPFfile) && isa(SDS, 'double') && isa(tissueType, 'char')
        if nargin < 7
            plotChannelIdx = [];
        end
    
        % Load the main data file with the intensity
        data = load(dataFile);
        Lambda = data.SD.Lambda;
        timeVector = data.t;
        intensityData = data.d;
        
        % Read the extinction coefficients from the CSV file and extract
        % the values for the Hb02 and the HHb for our further analysis.
        % This is done for both of the wavelengths of interest (both of the
        % lambdas.
        extinctionData = readtable(extinctionCoefficientsFile);
        extCoeffLambda1 = extinctionData{extinctionData.wavelength == Lambda(1), {'HbO2', 'HHb'}};
        extCoeffLambda2 = extinctionData{extinctionData.wavelength == Lambda(2), {'HbO2', 'HHb'}};
    
        % Read the DPF data specific for the tissue type that is inputted
        DPFData = readtable(DPFperTissueFile, 'Delimiter', '\t');
        DPF = DPFData{strcmp(DPFData.Tissue, tissueType), 'DPF'};
    
        % Read the relative DPF coefficients - this is specific for each
        % wavelength
        relDPFData = readtable(relDPFfile);
        relDPF1 = relDPFData{relDPFData.wavelength == Lambda(1), 'relDPFcoeff'};
        relDPF2 = relDPFData{relDPFData.wavelength == Lambda(2), 'relDPFcoeff'};
    
        % We then multiply the wavelength DPF (relative DPF) by the DPF
        % specific to the tissue, to get the DPF that we will use in our
        % further calculations
        totalDPF1 = DPF * relDPF1;
        totalDPF2 = DPF * relDPF2;
    
        % Split intensity data into two wavelengths - the first 20 columns
        % are the first wavelength and the last 20 columnds are the second
        % wavelength
        I1 = intensityData(:, 1:20);
        I2 = intensityData(:, 21:40);
    
        % In order to calculate the change in optical density we take the log
        % difference between the intensity at the beginning with the
        % intensity at every time point - we do this for each wavelength
        % again
        OD1 = log(I1);
        OD2 = log(I2);
        dOD1 = [-OD1 + OD1(1,:)];
        dOD2 = [-OD2 + OD2(1,:)];
    
        % The total pathlength is the total DPF multiplied by the length
        % it has to travel which is the SDS - the source detector
        % separation
        pathlength_1 = totalDPF1 * SDS;
        pathlength_2 = totalDPF2 * SDS;
    
        % We created a matrix that calculates for each wavelength and for
        % each pathlength, the extinction coefficient so we could then
        % multiply it by the intensity and get all of the values for each
        % of the wavelengths
        lambda1HbO = extCoeffLambda1(1);
        lambda1HbR = extCoeffLambda1(2);
        lambda2HbO = extCoeffLambda2(1);
        lambda2HbR = extCoeffLambda2(2);
        
        E = [lambda1HbO/pathlength_1, lambda1HbR/pathlength_1; ...
             lambda2HbO/pathlength_2, lambda2HbR/pathlength_2];
    
        % We take the inverse (as shown in the powerpoint)
        invE = inv(E);
    
        % We initiate a matrix to calculate for all channels at every type
        % point
        dHbO = zeros(size(dOD1));
        dHbR = zeros(size(dOD1));
    
        % For each of the wavelengths (represented by the dOD1 and dOD2
        % calculations that we did) we calculate the coefficients multiplied by  the
        % intensity to get the blood concentration at every point for every
        % column 
        for i = 1:20
            % We take apart the matrix and assign the values to the
            % relevant HbO part and HbR part
            dOD = [dOD1(:,i), dOD2(:,i)];
            dC = invE * dOD';
            dHbO(:, i) = dC(1,:);
            dHbR(:, i) = dC(2,:);
        end
    
        % Plot if required - only if we give channel indices to plot
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
    else
        disp("There are missing parameters or parameters of the wrong type")
    end
end

