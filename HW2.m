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
ProjDataRamp = fftshift(fft(ProjDataWeighted,[],2),2);
ProjDataRamp = ProjDataRamp .* rampFilter;
ProjDataRamp = ifftshift(ProjDataRamp,2);
ProjDataRamp = real(ifft(ProjDataRamp,[],2));
save('ProjDataRamp.mat', 'ProjDataRamp');

imshow(ProjDataRamp, []);
title('Filtered Projection Data');

ReconstructedPhantom = zeros(phantomLen);
Xcenter = [128,128]; %array center index of phantom
X = [0,0]; %(x,y) coords of phantom

%Source location and unit vectors for each angle
SrcLocs = (DistD/2 * [cos(Angle' ), sin(Angle' )]);
d1 = [-cos(Angle' ), -sin(Angle' )];
d2 = [-sin(Angle' ), cos(Angle' )];

deltaX = 2*(Radius)/phantomLen; %pixel width
deltaBeta = 2*pi/YL; %angular resolution

%for each pixel, compute offset from center
for X1 = 1:256
    X(1) = (Xcenter(1) - X1)*deltaX;
    for Y1 = 1:256
        X(2) = (Y1 - Xcenter(2))*deltaX;
        
        %For each angle, backproject data to pixels
        for AngIdx = 1:720
        
            %compute distance of point from source
            SrcLoc = SrcLocs(AngIdx,:);
            PtSrcDiffVector = SrcLoc - X;
            
            %compute projection onto detector
            u1 = -DistD * PtSrcDiffVector * d2(AngIdx,:)';
            u1 = u1 / (PtSrcDiffVector * d1(AngIdx,:)');

            D = norm(PtSrcDiffVector);
            
            %find closest detector index
            projIndex = round((u1/deltaU)+300);
            
            %Some points fall outside of detector. Do not backproject them.
            %Typically these are in the corners of the image. All phantom
            %pts found
            if((projIndex > 0) && (projIndex < YL + 1))
                %Weight and perform backprojection
                ProjValue = deltaBeta * ((DistD^2 + u1^2)*50/((D^2)*DistD)) * ProjDataRamp(AngIdx,projIndex);
                ReconstructedPhantom(257-Y1,257-X1) = ReconstructedPhantom(257-Y1,257-X1) + ProjValue;
            end
        end
    end
end
ReconstructedPhantom = 0.5 .* ReconstructedPhantom;
figure;
imshow(ReconstructedPhantom, []);
title('Reconstructed Phantom');
save('ReconstructedPhantom.mat', 'ReconstructedPhantom');
