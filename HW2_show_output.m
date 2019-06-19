load('ReconstructedPhantom.mat')
load('ProjDataRamp.mat')

subplot(1,2,1); imshow(ProjDataRamp,[]); title('Filtered Projection Data');
subplot(1,2,2); imshow(ReconstructedPhantom,[]); title('Reconstructed Phantom');
