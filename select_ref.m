function ref_num = select_ref(img,repnum,echon,threshold)
% Remove images with artifacts.
blur = blur_measure(img,repnum,echon,threshold);
normal = setdiff(1:10,blur);
while isempty(normal)
    threshold = threshold - 0.01;
    blur = blur_measure(img,repnum,echon,threshold);
    normal = setdiff(1:10,blur);
end
S = zeros(length(normal),length(normal));
for i = 1:length(normal)
    for j = 1:length(normal)
    dis = img(:,:,normal(i),:) - img(:,:,normal(j),:);
    S(i,j) = norm(dis(:));
    end
end

S = sum(S,2);
mins = find(S == min(S));
ref_num = normal(mins(1));
end

