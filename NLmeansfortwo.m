function [dist,signal] = NLmeansfortwo(ref,test,block_size,search_window,TEnum,block_num,maxm)
ref = ref / maxm;
testc = test / maxm;
% The core of Non-local means algorithm.
% dist: The weighted Euclidean distance
% signal: Shifted images
[len, width,~] = size(ref);
f = (block_size - 1) / 2;
t = (search_window - 1) / 2;

kernel = make_kernel(f);  %Gaussian kernel used to weighted distance of two blocks
kernel = kernel / sum(sum(kernel));
kernel = repmat(kernel,1,1,TEnum);
signal = zeros(len,width,1,block_num,TEnum);
dist = zeros(len,width,block_num);
testc = padarray(testc,[f f],'symmetric');
test = padarray(test,[f f],'symmetric');
ref = padarray(ref,[f f],'symmetric');
for i = 1:len
    for j = 1:width
        i1 = i + f;
        j1 = j + f; 
        block1 = ref(i1 - f:i1 + f,j1 - f:j1 + f,:);
        mmin = max(i1-t,f+1);
        mmax = min(i1+t,len+f);
        nmin = max(j1-t,f+1);
        nmax = min(j1+t,width+f);
        len1 = length(mmin:mmax);
        width1 = length(nmin:nmax);
        dist_record = zeros(len1,width1);
        for n = nmin:nmax
            for m = mmin:mmax
                a = m - f; b = m + f;
                x = n - f; y = n + f;
                d = block1 - testc(a:b,x:y,:);
                d = real(d).^2 + imag(d).^2;
                d = kernel .* d;
                a = m - mmin + 1; b = n - nmin + 1;
                dist_record(a, b) = sum(sum(sum(d)));% + sqrt((m - i1).^2 + (n - j1)^2);
            end
        end
        [b,a] = sort(dist_record(:));
        a = a(1:block_num); b = b(1:block_num);
        [x, y] = ind2sub([len1,width1],a);
        x = x + mmin - 1; 
        y = y + nmin - 1;   
        a = sub2ind([len + 2 * f,width + 2 * f],x,y);
        dist(i,j,:) = b;
        for k = 1:TEnum
            tmp = squeeze(test(:,:,k));
            signal(i,j,1,:,k) = tmp(a);
        end
    end
end


function [kernel] = make_kernel(f)
%% Make Gaussian kernel
kernel=zeros(2*f+1,2*f+1);
for d=1:f
    value= 1 / (2*d+1)^2 ;
    %value = 1 / (1.5*d+1)^2;
    for j=-d:d
        for i=-d:d
            kernel(f+1-i,f+1-j)= kernel(f+1-i,f+1-j) + value ;
        end
    end
end
kernel = kernel ./ f;
