function u = location(img)
[m,~,echon] = size(img);
img = abs(img);
S = zeros(1,echon);
for i = 1:echon
  img(:,:,i) = img(:,:,i) / max(max(squeeze(img(:,:,i))));
end
e = zeros(m,echon);
i = 1;
for j = 1:m
  e(j,i) = entropy(squeeze(double(img(j,:,i))));
end
[y,x] = findpeaks(squeeze(e(:,i)));
yy = y(2:end) - y(1:end-1);
xmax = find(yy == max(yy)) + 1;
u = x(xmax);
% img_u = img(u - 10:u,:,:);
% for i = 1:echon
%   img_u(:,:,i) = (img_u(:,:,i) / max(max(squeeze(img_u(:,:,i))))) * 100;
%   S(i) = estimate_sharpness(img_u(:,:,i));
% end
% display(S / min(S));
% figure; imagescn(squeeze(abs(img(:,:,:))),[],[2 3],[],2);
% measure = min(S);



