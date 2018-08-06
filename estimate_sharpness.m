function sharpness = estimate_sharpness(img)

[Gx,Gy] = gradient(img,1,0.1);
S = sqrt(Gx.*Gx + Gy.* Gy);
sharpness = sum(S(:))./(numel(Gx));
