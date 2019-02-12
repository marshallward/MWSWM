function [x, y, t, u, v, h, Lx] = readMyData(filename)

fid = fopen(filename,'r','ieee-le');

% Lx - physical x box length
size = fread(fid, 1, 'int32');
Lx = fread(fid, 1, 'double');
size = fread(fid, 1, 'int32');

% Ly - physical y box length
size = fread(fid, 1, 'int32');
Ly = fread(fid, 1, 'double');
size = fread(fid, 1, 'int32');

% dx - physical distance between x nodes
size = fread(fid, 1, 'int32');
dx = fread(fid, 1, 'double');
size = fread(fid, 1, 'int32');

% dy - physical distance between y nodes
size = fread(fid, 1, 'int32');
dy = fread(fid, 1, 'double');
size = fread(fid, 1, 'int32');

% dt - physical time between timesteps
size = fread(fid, 1, 'int32');
dt = fread(fid, 1, 'double');
size = fread(fid, 1, 'int32');

% nM - number of x nodes
size = fread(fid, 1, 'int32');
nM = fread(fid, 1, 'int32');
size = fread(fid, 1, 'int32');

% nN - number of y nodes
size = fread(fid, 1, 'int32');
nN = fread(fid, 1, 'int32');
size = fread(fid, 1, 'int32');

% nT - number of timesteps
size = fread(fid, 1, 'int32');
nT = fread(fid, 1, 'int32');
size = fread(fid, 1, 'int32');

% mM - maximum x node
size = fread(fid, 1, 'int32');
mM = fread(fid, 1, 'int32');
size = fread(fid, 1, 'int32');

% mN - maximum y node
size = fread(fid, 1, 'int32');
mN = fread(fid, 1, 'int32');
size = fread(fid, 1, 'int32');

% epsilon - Rossby number (nonlinear scale factor)
size = fread(fid, 1, 'int32');
epsilon = fread(fid, 1, 'double');
size = fread(fid, 1, 'int32');

% dataRate - Number of steps per recorded step
size = fread(fid, 1, 'int32');
dataRate = fread(fid, 1, 'int32');
size = fread(fid, 1, 'int32');

% buffRate - Number of stored steps per write
size = fread(fid, 1, 'int32');
buffSize = fread(fid, 1, 'int32');
size = fread(fid, 1, 'int32');

% x - physical x grid
size = fread(fid, 1, 'int32');
x = fread(fid, nM, 'double');
size = fread(fid, 1, 'int32');

% y - physical y grid
size = fread(fid, 1, 'int32');
y = fread(fid, nN, 'double');
size = fread(fid, 1, 'int32');

% t - physical y grid
size = fread(fid, 1, 'int32');
t = fread(fid, (nT/dataRate)+1, 'double');
size = fread(fid, 1, 'int32');

u = zeros(nM, nN, (nT/dataRate)+1);
v = zeros(nM, nN, (nT/dataRate)+1);
h = zeros(nM, nN, (nT/dataRate)+1);

% Get the initial states
size = fread(fid, 1, 'int32');
u(:,:,1) = reshape(fread(fid, nM*nN, 'double'), nM, nN);
size = fread(fid, 1, 'int32');

size = fread(fid, 1, 'int32');
v(:,:,1) = reshape(fread(fid, nM*nN, 'double'), nM, nN);
size = fread(fid, 1, 'int32');

size = fread(fid, 1, 'int32');
h(:,:,1) = reshape(fread(fid, nM*nN, 'double'), nM, nN);
size = fread(fid, 1, 'int32');

% Grabbin da dataBuffer
for i = 1:(nT/dataRate/buffSize)
    size = fread(fid, 1, 'int32');
    data = reshape(fread(fid, nM*nN*buffSize*3, 'double'), nM, nN, buffSize, 3);
    size = fread(fid, 1, 'int32');

    for j = 1:buffSize
        u(:,:,i*j+1) = squeeze(data(:,:,j,1));
        v(:,:,i*j+1) = squeeze(data(:,:,j,2));
        h(:,:,i*j+1) = squeeze(data(:,:,j,3));
    end
end

fclose(fid);