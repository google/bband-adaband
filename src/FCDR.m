% Copyright 2021 Google LLC
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%      http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

function [ fcc_map_1, fcc_map_2, fcc_map_3 ] = FCDR(img, th1, th2, th3, th4, th5)
%% implemented fcdr algorithm
% ref:
% Huang, Qin, et al. "Understanding and removal of false contour in hevc 
% compressed images." IEEE Transactions on Circuits and Systems for Video
% Technology 28.2 (2016): 378-391.
%
% Input:
%
% img: 2D color (rgb) or intensity image
% th1-th5: thresholds defined in the paper. Already tuned. Not recommended to
% change.
%
% Output:
%
% fcc_map: false contour candidates map
% fcc_map_3: final false contour map
%
% Required toolboxes: None
%

%% set default params
if nargin<=1
  % Step1
  th1 = 1; % Smooth are thresh,               bigger, tighter
  % Step2
  th2 = 3; % textures thresh, max,            smaller, tighter
  th3 = 7; % textures thresh, sum,            smaller, tighter
  % Step3
  th4 = 10; % monotonicity thresh, max,        bigger, tighter
  % fixed
  th5 = 6; % monotonicity thresh, min,        smaller, tighter
end

% convert to gray image if input is not
[~,~,z] = size(img);
if z==3
  img = rgb2gray(img);
else
  assert(z==1, 'unknown image type. Input image should be grayscale or rgb');
end

img = double(img);

[h, w] = size(img);
% initialize fcc_map
fcc_map_0 = ones(size(img));

% replicate padding for gradient comp
img_pad = padarray(img, [1 1], 'replicate', 'both');

% calculate gradient among 8 directions

% 8 dirs: [tl, top, tr, right, br, bottom, bl, left]
%         [1,   2,   3,  4,    5,     6,   7,    8]
grad_dirs = [-1, -1; -1, 0; -1, 1; 0, 1; 1, 1; 1, 0; 1, -1; 0, -1];
% initialize gradient map
grad_map = zeros([size(img), size(grad_dirs, 1)]);

for i=1:h % row
  for j=1:w % col
    ni = i+1; % updated idx in padded image
    nj = j+1; % updated idx in padded image
    for dir=1:size(grad_dirs, 1)
      % calculate grad for each direction
      grad_map(i,j,dir) = abs(img_pad(ni,nj) - ...
      img_pad(ni+grad_dirs(dir,1), nj+grad_dirs(dir,2)));
    end
  end
end
% sum of gradients in opposite directions
grad_sum_15 = grad_map(:,:,1) + grad_map(:,:,5);
grad_sum_26 = grad_map(:,:,2) + grad_map(:,:,6);
grad_sum_37 = grad_map(:,:,3) + grad_map(:,:,7);
grad_sum_48 = grad_map(:,:,4) + grad_map(:,:,8);
% concat into one nd matrix
grad_sum_mat = zeros([size(img), size(grad_dirs,1)/2]);
grad_sum_mat(:,:,1) = grad_sum_15;
grad_sum_mat(:,:,2) = grad_sum_26;
grad_sum_mat(:,:,3) = grad_sum_37;
grad_sum_mat(:,:,4) = grad_sum_48;
% calculate grad mag and direction

grad_max_mat = max(grad_sum_mat, [], 3);
grad_mag_mat = sum(grad_sum_mat, 3);

% Step1: remove smooth regions
fcc_map_1 = fcc_map_0 & (grad_sum_15 >= th1) | (grad_sum_26 >= th1) |...
    (grad_sum_37 >= th1) | (grad_sum_48 >= th1);

% Step2: exclude textures and sharp edges
fcc_map_2 = fcc_map_1 & (grad_max_mat < th2) ...
    & (grad_mag_mat < th3);

% Step3: exclude areas without monotonicity
% slightly modified:
% given a edge, the 5x5 neighors, count number of similar-gradients
% pixels along each direction
fcc_map_3 = fcc_map_2;
blksize_row = 9; % processing in the neighboring block centered at this pixel
blksize_col = 9;
padding_row = floor(blksize_row/2);
padding_col = floor(blksize_col/2);
grad_sum_mat_pad = padarray(grad_sum_mat, [padding_row padding_col], ...
    'replicate', 'both'); % padding image for corner pixels
% find edges found in fcc_map_2
[ei,ej] = find(fcc_map_2);
% store # similar pixels in principal direction, for debug
sim_cnt_max_mat = zeros(1,length(ei));
% store # similar pixels in norm direction, for debug
sim_cnt_nor_mat = zeros(1,length(ei));
for idx = 1:length(ei) % only traverse the edges in fcc_map_2
  ii = ei(idx);
  jj = ej(idx); % get edge coord
  % extract block centered at (ii,jj)
  grad_sum_blk = grad_sum_mat_pad(ii:ii+blksize_row-1,jj:jj+blksize_col-1, :);
  sim_th = 2;
  grad_sim_blk = sum(abs(grad_sum_blk - grad_sum_blk(ceil(blksize_row/2),...
      ceil(blksize_col/2),:)), 3) < sim_th;
  sim_cnt_ver = sum(sum(grad_sim_blk(:,floor(blksize_col/2):ceil(blksize_col/2)+1)));
  sim_cnt_hon = sum(sum(grad_sim_blk(floor(blksize_col/2):ceil(blksize_row/2)+1,:)));
  sim_cnt_dia = sum(diag(grad_sim_blk)) + sum(diag(grad_sim_blk(2:end,1:end-1))) + ...
      sum(diag(grad_sim_blk(1:end-1,2:end))) ;
  grad_sum_blk_fl = fliplr(grad_sim_blk);
  sim_cnt_ant = sum(diag(grad_sum_blk_fl)) + sum(diag(grad_sum_blk_fl(2:end,1:end-1))) + ...
      sum(diag(grad_sum_blk_fl(1:end-1,2:end))) ;
  sim_cnt_max = max([sim_cnt_ver, sim_cnt_hon, ...
      sim_cnt_dia, sim_cnt_ant]);
  sim_cnt_max_mat(idx) = sim_cnt_max;
  sim_cnt_min = min([sim_cnt_ver, sim_cnt_hon, ...
      sim_cnt_dia, sim_cnt_ant]);
  if sim_cnt_max == sim_cnt_ver
    sim_cnt_n = sim_cnt_hon;
  elseif sim_cnt_max == sim_cnt_hon
    sim_cnt_n = sim_cnt_ver;
  elseif sim_cnt_max == sim_cnt_dia
    sim_cnt_n = sim_cnt_ant;
  else
    sim_cnt_n = sim_cnt_dia;
  end
  sim_cnt_nor_mat(idx) = sim_cnt_n;
  if sim_cnt_max < th4 || sim_cnt_min > th5
    fcc_map_3(ii,jj) = 0;
  end
end
end
