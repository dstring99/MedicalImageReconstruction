clear;
load('Res_Proj_Equi_Dist_Original_Phantom');

%Recontruct 256x256 phantom
phantomLen = 256;

%Compute cosine pre-weighting factor
deltaU = DecLength/YL;
u = deltaU * (-YL/2+1:YL/2);
weight = sqrt((DistD*DistD) + (u.*u));
weight = DistD ./ weight;
weight = repmat(weight, [720,1]);

%Weight projection data
ProjDataWeighted = ProjData .* weight;

%Create ramp filer in frequency domain
rampFilter = abs(linspace(-1, 1, 600));
rampFilter = repmat(rampFilter, [720,1]);

%Filter the weighted data
ProjDataFft = fft(ProjDataWeighted,[],2);
ProjDataRamp = fftshift(ProjDataFft,2);

figure(1);
subplot(1,3,1);
plot(rampFilter(90,:));
title('Ramp Filter');
subplot(1,3,2);
plot(abs(ProjDataFft(180,:)));
title('Projection Data FFT');
subplot(1,3,3);
plot(abs(ProjDataRamp(180,:)));
title('Projection Data FFT Shifted');

ProjDataRamp = ProjDataRamp .* rampFilter;
ProjDataRamp = ifftshift(ProjDataRamp,2);
ProjDataRamp = real(ifft(ProjDataRamp .* rampFilter,[],2));
save('ProjDataRamp.mat', 'ProjDataRamp');

figure(2);
subplot(1,2,1);
imshow(ProjData, []);
title('Projection Data');
subplot(1,2,2);
imshow(ProjDataRamp, []);
title('Filtered Projection Data');