name = 'PLP_P02';
nx = 4; %How many columns there are
ny = 96; %How many rows per column
dx = 16*10^-6*(10^-3); %The spacing in x direction in mm
dy = 20*10^-6*(10^-3); %The spacing in y direction in mm
h = 0.1; %The span of the CSD in z in mm
ycoords_mm = ycoords.*10^-3; %Convert y coordinates into mm
xcoords_mm = xcoords.*10^-3; %Convert x coordinates into mm
ix_chans = [1:n_func_chan]';
dsp(:,:,1) = [ix_chans,xcoords_mm]; % x-displacement
dsp(:,:,2) = [ix_chans,ycoords_mm]; % y-displacement
initlin2d(name,nx,ny,dx,dy,h,dsp);