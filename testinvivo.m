cBasepath = '/Users/huiwenluo/Desktop/data_cse2d_huiwen';
path = '/Users/huiwenluo/Desktop/data_cse2d_huiwen/added subjects';
file = dir(path);


for i = 9
    pn = file(i+3).name;
    pathn = dir([path,'/',pn]);
    data = load([path,'/',pn,'/',pathn(3).name]);
    imscc = data.imscc;
    TE = data.TE;
    Size = size(imscc);
    len = Size(1); width = Size(2);
    Slice_num = Size(4);  Rnum = Size(3); TEnum = Size(6);
    %select the references using edge detection
    ref_num = ones(length(file)-2,Slice_num);
    for j = 1:Slice_num

    display(['slice:',num2str(j)]);
    ref_num(i,j) = selcect_ref(squeeze(imscc(:,:,:,j,1,:)),Rnum,TEnum,0.25);
    end
    if mod(len,2) ~= 0
        len = len - 1;
        imscc = imscc(2:end,:,:,:,:,:);
    end
    if mod(width,2)~=0
        width =width - 1;
        imscc = imscc(:,2:end,:,:,:,:);
    end
    output_nlm = zeros(len,width,Slice_num,TEnum);
    r2starmap_nlm = zeros(len,width,Slice_num);
    fatfraction_nlm = zeros(len,width,Slice_num);
    r2starmap_avg = zeros(len,width,Slice_num);
    fatfraction_avg = zeros(len,width,Slice_num);
    output_avg = squeeze(mean(imscc,3));

    for j = 1:Slice_num 
%         output_ref(:,:,j,:) = squeeze(imscc(:,:,ref_num(i,j),1,:));
%     [r2starmap_avg(:,:,j),fatfraction_avg(:,:,j)] = processFatWater((reshape(output_avg(:,:,j,:),len,width,1,1,1,TEnum)) ,TE(1),TE(2) - TE(1),TEnum);
    [output_nlm(:,:,j,:),w,s,output_ref(:,:,j,:)] = Process(ref_num(i,j),j,imscc,TE);%ref_num(i),i,imscc,TE);
    end
    save([path,'/',pn,'/result.mat'],'output_nlm','output_avg','output_ref','TE');
  
end





