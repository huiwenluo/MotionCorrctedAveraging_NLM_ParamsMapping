function blur = blur_measure(imscc,repnum,echon,threshold)

% u = blur_measure(squeeze(imscc(:,:,i,6,1,:)));
imscc = squeeze(imscc);
for i = 1:repnum
    m(i) = location(squeeze(imscc(:,:,i,:)));
end

u = max(m);
img = abs(imscc(u - 15:u + 2,:,:,:));
S = zeros(1,echon);
blur = [];
for k = 1:repnumc
    img_u = squeeze(img(:,:,k,:));
    for i = 1:echon
        img_u(:,:,i) = (img_u(:,:,i) / max(max(squeeze(img_u(:,:,i)))));
        tmp = 1 - detectedges(((img_u(:,:,i))));
        S(i) = estimate_sharpness(tmp);
        if S(i) < threshold
            blur = [blur,k];
            break;
        end
    end
%     display(S);
%     imagescn(squeeze(abs(imscc(:,:,k,:))),[],[2 3],[],[]);
%     close;
end
