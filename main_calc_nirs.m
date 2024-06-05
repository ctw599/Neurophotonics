main_data_file = "C:\Users\user\Downloads\HW1_data\FN_032_V1_Postdose1_Nback.mat";
SDS_cm = 3;
tissue_type = 'adult_forearm';
extinction_coefficients_file = "C:\Users\user\Downloads\HW1_data\ExtinctionCoefficientsData.csv";
DPF_per_tissue_file = "C:\Users\user\Downloads\HW1_data\DPFperTissue.txt";
rel_DPF_file = "C:\Users\user\Downloads\HW1_data\RelativeDPFCoefficients.csv";
plot_channel_idx = [1, 2];
[dHbR, dHbO, fig] = CalcNIRS(main_data_file, SDS_cm, tissue_type, extinction_coefficients_file, DPF_per_tissue_file, rel_DPF_file, plot_channel_idx);