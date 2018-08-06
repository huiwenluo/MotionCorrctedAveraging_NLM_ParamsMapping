function [output,w,s] = nlmslice(ref,imscc,i,block_size,search_window,block_num,slice_window,ref_num,ss)
Size = size(imscc);
len = Size(1); width = Size(2);
Slice_num = Size(4);  Rnum = Size(3); TEnum = Size(6);
%% just remove the reference image to perform NLM better.
if i <= (slice_window - 1)/2
    order = 1:slice_window;
    order(order == i) = [];
    group = reshape(imscc(:,:,:,order,1,:),len,width,Rnum * (slice_window - 1),TEnum);
    if ref_num == 1
        order = (2:Rnum);
    else
    order = [1:ref_num - 1, ref_num + 1 : Rnum];
    end
    group = cat(3,group,reshape(imscc(:,:,order,i,1,:),len,width,Rnum - 1,TEnum));
elseif i >= Slice_num - (slice_window - 1)/2
    order = Slice_num - (slice_window - 1 ):Slice_num;
    order(order == i) = [];
    group = reshape(imscc(:,:,:,order,1,:),len,width,Rnum * (slice_window - 1),TEnum);
    order = [1:ref_num - 1, ref_num + 1 : Rnum];
    group = cat(3,group,reshape(imscc(:,:,order,i,1,:),len,width,Rnum - 1,TEnum));
else
    order = i - (slice_window - 1) / 2 : i + (slice_window - 1) / 2;
    order(order == i) = [];
    group = reshape(imscc(:,:,:,order,1,:),len,width,Rnum * (slice_window - 1),TEnum);
    order = [1:ref_num - 1, ref_num + 1 : Rnum];
    group = cat(3,group,reshape(imscc(:,:,order,i,1,:),len,width,Rnum - 1,TEnum));
end


%%  Compele algorithm

sigma = 0.0001;
dist = squeeze(zeros(len,width,(Rnum * slice_window - 1) , block_num));
signal = zeros(len,width,(Rnum * slice_window - 1), TEnum);
maxm = max(abs(imscc(:)));
parfor k = 1:(Rnum * slice_window - 1)
    [dist(:,:,k),signal(:,:,k,:)] = NLmeansfortwo(ref,squeeze(group(:,:,k,:)),block_size,search_window,TEnum,block_num,maxm);
end

weight = exp(-dist/(sigma * ss^2));
%weight = exp(dist/(sigma * ss^2));
% weight = dist;
wref = ones(len,width);  % Set weights of reference image to one
weight = cat(3,weight,wref);
weight = weight ./ repmat(sum(weight,3),1,1,(Rnum * slice_window - 1) * block_num + 1);
w = weight;
weight = repmat(weight,1,1,1,TEnum);
signal = cat(3,signal,reshape(ref,len,width,1,TEnum));
s = signal;
output = sum(weight .* signal,3);
