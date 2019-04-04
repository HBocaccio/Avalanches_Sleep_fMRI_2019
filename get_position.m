function output = get_position(i,j,k,xmax,ymax,zmax)
% Get single-index target from 3D position
output = ymax*zmax*(i-1) + zmax*(j-1) + k;
end   