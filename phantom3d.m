function [p,ellipse]=phantom3d(slice_num,respiZ,respiO,n,Amp,varargin)
%   Generate the slice_num-th slice of abdomen phantom. 
%   slice_num: the sequence number of slice
%   respiZ: variation of the length of axis of the ellipsoid to change the shape
%   of organs.(Anterior/Posterior)
%   respiO: the same function as respiZ. (Superior/Inferior)
%   Parameters of the ellipsoid has ten columns,
%   with each column containing a different parameter for the ellipsoids:
%
%     Column 1:  A      the additive intensity value of the ellipsoid
%     Column 2:  a      the length of the x semi-axis of the ellipsoid
%     Column 3:  b      the length of the y semi-axis of the ellipsoid
%     Column 4:  c      the length of the z semi-axis of the ellipsoid
%     Column 5:  x0     the x-coordinate of the center of the ellipsoid
%     Column 6:  y0     the y-coordinate of the center of the ellipsoid
%     Column 7:  z0     the z-coordinate of the center of the ellipsoid
%     Column 8:  phi    phi Euler angle (in degrees) (rotation about z-axis)
%     Column 9:  theta  theta Euler angle (in degrees) (rotation about x-axis)
%     Column 10: psi    psi Euler angle (in degrees) (rotation about z-axis)
%

ellipse = parse_inputs(respiZ,respiO);
p = zeros(n,n);
rng =  ( (0:n-1)-(n-1)/2 ) / ((n-1)/2);

[x,y,z] = meshgrid(rng,rng,rng);
x = x(:,:,slice_num);
y = y(:,:,slice_num);
z = z(:,:,slice_num);
coord = [flatten(x); flatten(y); flatten(z)];
p = flatten(p);

for k = 1:size(ellipse,1)
    A = Amp(k);      % Amplitude change for this ellipsoid
    asq = ellipse(k,2)^2;        % a^2
    bsq = ellipse(k,3)^2;        % b^2
    csq = ellipse(k,4)^2;        % c^2
    x0 = ellipse(k,5);           % x offset
    y0 = ellipse(k,6);           % y offset
    z0 = ellipse(k,7);           % z offset
    phi = ellipse(k,8)*pi/180;   % first Euler angle in radians
    theta = ellipse(k,9)*pi/180; % second Euler angle in radians
    psi = ellipse(k,10)*pi/180;  % third Euler angle in radians
    
    cphi = cos(phi);
    sphi = sin(phi);
    ctheta = cos(theta);
    stheta = sin(theta);
    cpsi = cos(psi);
    spsi = sin(psi);
    
    % Euler rotation matrix
    alpha = [cpsi*cphi-ctheta*sphi*spsi   cpsi*sphi+ctheta*cphi*spsi  spsi*stheta;
        -spsi*cphi-ctheta*sphi*cpsi  -spsi*sphi+ctheta*cphi*cpsi cpsi*stheta;
        stheta*sphi                  -stheta*cphi                ctheta];
    
    % rotated ellipsoid coordinates
    coordp = (alpha*coord);
    coordc = coordp(1:2,:);
    idx = find((coordc(1,:)-x0).^2./asq + (coordc(2,:)-y0).^2./bsq  <= 1 -(coordp(3,:)-z0).^2./csq);
    if k == 3
        flag = idx;
    elseif k == 4
        idx = setdiff(idx,intersect(idx,flag));
    end
    p(idx) = p(idx) + A;
    
end

p = reshape(p,[n n]);

return;


function out = flatten(in)

out = reshape(in,[1 prod(size(in))]);

return;


function e = parse_inputs(respiZ,respiO)
%  e is the m-by-10 array which defines ellipsoids
%  n is the size of the phantom brain image
% The default size
e = modified_shepp_logan(respiZ,respiO);
% ellipse is not yet defined
if isempty(e)
    e = modified_shepp_logan(respiZ,respiO);
end

return;


function e = modified_shepp_logan(respiZ,respiO)
%     Column 1:  A      the additive intensity value of the ellipsoid
%     Column 2:  a      the length of the x semi-axis of the ellipsoid
%     Column 3:  b      the length of the y semi-axis of the ellipsoid
%     Column 4:  c      the length of the z semi-axis of the ellipsoid
%     Column 5:  x0     the x-coordinate of the center of the ellipsoid
%     Column 6:  y0     the y-coordinate of the center of the ellipsoid
%     Column 7:  z0     the z-coordinate of the center of the ellipsoid
%     Column 8:  phi    phi Euler angle (in degrees) (rotation about z-axis)
%     Column 9:  theta  theta Euler angle (in degrees) (rotation about x-axis)
%     Column 10: psi    psi Euler angle (in degrees) (rotation about z-axis)
%
%         A      a     b     c     x0      y0      z0    phi  theta    psi
%        -----------------------------------------------------------------
e =    [  1   .920 .6900  inf      0       0       0      0      0      0
    -.8   .900 .6700  inf      0       0       0      0      0      0
    1  .3000  .510  .620    .30   0.184      0     -50   -30     10
    1  .3000  .510  .620   -.30    .184      -0.2      50    30     10
    .1  .100   .100  inf      0      .35   -.15      0      0      0 ];
e(3:4,7) = e(3:4,7) + respiZ;
e(3:4,2) = e(3:4,2) + respiO;

return;

