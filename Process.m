function [output,w,s,ref] = Process(ref_num,i,imscc,TE)%1,dTE)%,R2starmap,fatfraction2] = Process(ref_num,i,imscc,TE)%1,dTE)
Size = size(imscc);
len = Size(1); width = Size(2);
Slice_num = 11;  Rnum = Size(3); TEnum = Size(6);
if mod(len,2) ~= 0
    len = len - 1;
    imscc = imscc(2:end,:,:,:,:,:);
end
if mod(width,2) ~= 0
    width = width - 1;
    imscc = imscc(:,2:end,:,:,:,:);
end
% weight threshold?
ref = squeeze(imscc(:,:,ref_num,i,1,:));
%%  h
ss = 10;
slice_window = 3;
block_size = 11;
search_window = 3;
block_num = 1;
[output,w,s] = nlmslice(ref,imscc,i,block_size,search_window,block_num,slice_window,ref_num,ss);
% [outParams2,fatfraction2] = processFatWater((reshape(output(:,:,:,:),len,width,1,1,1,TEnum)) ,TE(1),TE(2) - TE(1),TEnum);
%  R2starmap = outParams2.r2starmap;
