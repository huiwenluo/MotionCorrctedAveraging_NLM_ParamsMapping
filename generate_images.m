function image = generate_images(sample_f,sample_w)
%% Generate real images without noise by sample of fat-fraction map and water-fraction map.

%% Add to matlab path
BASEPATH = '~/Desktop/code_cse2d_huiwen/012_graphcut3D/';
addpath([BASEPATH 'common/']);
addpath([BASEPATH 'graphcut/']);
addpath([BASEPATH 'descent/']);
addpath([BASEPATH 'mixed_fitting/']);
addpath([BASEPATH 'create_synthetic/']);
addpath([BASEPATH 'matlab_bgl/']);

sx = 144;
sy = 144;
[sx1,sy1] = size(sample_f);
f_sample_f = fftshift(fftshift(fft(fft(ifftshift(ifftshift(sample_f,1),2),[],1),[],2),1),2)/sqrt(sx1*sy1);
f_sample_w = fftshift(fftshift(fft(fft(ifftshift(ifftshift(sample_w,1),2),[],1),[],2),1),2)/sqrt(sx1*sy1);
pick_x = (floor(sx1/2)+1- floor(sx/2)):(floor(sx1/2)+ floor(sx/2));
pick_y = (floor(sy1/2)+1- floor(sy/2)):(floor(sy1/2)+ floor(sy/2));
f_sample_f_2 = f_sample_f(pick_x,pick_y);
sample_f_2 = fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(f_sample_f_2,1),2),[],1),[],2),1),2)*sqrt(sx*sy);
f_sample_w_2 = f_sample_w(pick_x,pick_y);
sample_w_2 = fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(f_sample_w_2,1),2),[],1),[],2),1),2)*sqrt(sx*sy);





%% 
%trueParams.species(1).amps =  sample_w_2;% Water
%trueParams.species(2).amps =  sample_f_2;% Fat

trueParams.species(1).amps =  sample_w;% Water
trueParams.species(2).amps =  sample_f;% Fat


trueFF = computeFF(trueParams);

truer2star = 0*trueFF;
truer2star(trueFF>75) = 20;
truer2star(trueFF<75 & trueFF > 25) = 120;
truer2star(trueFF<25) = 70;
truer2star(abs(sample_w)+abs(sample_f)<15) = 0;

trueParams.fieldmap = zeros(sx1,sy1);
trueParams.r2starmap = truer2star;


%% Set data parameters
N = 8; TEinit = 1.2e-3; dTE = 1.4e-3;
TE = TEinit + [0:N-1]*dTE;
imDataParams0.TE = TE;
imDataParams0.FieldStrength = 3.0;
imDataParams0.PrecessionIsClockwise = -1;

%% Set recon parameters
% General parameters
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 0;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];


%% Simulate data
imDataParams1 = createSynthetic_imageSpace( imDataParams0, algoParams, trueParams );



fim = fftshift(fftshift(fft(fft(ifftshift(ifftshift(imDataParams1.images,1),2),[],1),[],2),1),2)/sqrt(sx1*sy1);
pick_x = (floor(sx1/2)+1- floor(sx/2)):(floor(sx1/2)+ floor(sx/2));
pick_y = (floor(sy1/2)+1- floor(sy/2)):(floor(sy1/2)+ floor(sy/2));
fim_2 = fim(pick_x,pick_y,:,:,:);
imDataParams.images = fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(fim_2,1),2),[],1),[],2),1),2)*sqrt(sx*sy);
imDataParams.images = 100*imDataParams.images/max(abs(imDataParams.images(:)));
imDataParams0 = imDataParams;

%%
% SNR = 10;
% noise_std = max(abs(imDataParams.images(:)))/SNR;
% imDataParams.images = imDataParams.images + noise_std*complex(randn(size(imDataParams.images)),randn(size(imDataParams.images)));

image = imDataParams.images;