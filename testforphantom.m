%% 
n = 512; 
respirFreq = 5/60;
respiZampl = 0.08;
respiOampl = 0.05;
A_w = [0,0,90,90,50];
A_f = [0,100,-90,-90,-50];
Rep_num = 10;
%%   generate real images without noises
t = 0.2:1:30;
respiZ = respiZampl * cos(2*pi * respirFreq * t).^(2*3);
respiO = respiOampl * cos(2*pi * respirFreq * t).^(2*3);
k = 1;
z_loc = 200:30:265;
thickness = 7;
imscc = zeros(144,144,Rep_num,length(z_loc),8);
img_f = zeros(n,n,thickness);
img_w = zeros(n,n,thickness);
tmp = -(thickness - 1)/2 :(thickness - 1)/2;
for i = 1:length(z_loc)
  for j = 1:Rep_num
        for h = 1:thickness
           img_w(:,:,h) = phantom3d(z_loc(i) + tmp(h),respiZ(k),respiO(k),n,A_w);
           img_f(:,:,h) = phantom3d(z_loc(i) + tmp(h),respiZ(k),respiO(k),n,A_f);
        end
         img = generate_images(sum(img_f,3),sum(img_w,3)); 
        imscc(:,:,j,i,:) = generate_images(sum(img_f,3),sum(img_w,3)); 
        k = k + 1;
    end
end
%%  Statistical process of Non-local means
imscc = reshape(imscc,144,144,10,3,1,8);
SNR = 10;
N = 8; TEinit = 1.2e-3; dTE = 1.4e-3;
TE = TEinit + [0:N-1]*dTE;
N = 50;    
r2starmap_nlm = zeros(144,144,N);
fatfraction_nlm = zeros(144,144,N);
output_nlm = zeros(144,144,8,N);
s = zeros(144,144,30,8,N);
w = zeros(144,144,30,N);
image = zeros(144,144,10,3,1,8,N);
for i = 1:N
    noise_std = max(abs(imscc(:)))/SNR;
    image(:,:,:,:,:,:,i) = imscc + noise_std*complex(randn(size(imscc)),randn(size(imscc)));
    tmp = respiZ(10:20);
    tmp = find(tmp == max(tmp));
    ref_num = tmp(1);
%     %[output_nlm(:,:,:,i),w(:,:,:,i),s(:,:,:,:,i),ref] = Process(ref_num,2,image(:,:,:,:,:,:,i),TE); 
     [output_nlm(:,:,:,i),w(:,:,:,i),s(:,:,:,:,i),ref,r2starmap_nlm(:,:,i),fatfraction_nlm(:,:,i)] = Process(ref_num,2,image(:,:,:,:,:,:,i),TE); 
end



