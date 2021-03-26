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

function [ band_vis_map, band_edge_list, band_score, grad_mag, grad_dir] = ...
    BBAD_I ( im_y, imguidedfilt_ws, thr1, thr2 )
%% Blind Banding Artifact Detector for image (BBAD_I)
%
% [ band_vis_map, band_edge_list, band_score, grad_mag, grad_dir] = ...
%     BBAD_I ( im_y, imguidedfilt_ws, thr1, thr2 )
% Detect the banding edge visibility map, banding edge map and banding score in
% an image
%
% Input Arguments
%
% im_y: input 1-channel image, could be grayscale, or Y-plane of YUV
% imguidedfilt_ws: window size for image guided pre-filtering. The larger window
% size, the smoother detected edge map obtained. Default value 1 (no smoothing).
% 1 is recommended to perform on VP9-transcoded frame, 5-9 is recommended to
% perform on RAW frame or H264/HEVC transcoded frame.
%
% thr1: threshold for extracting flat regions. Default value 2. Not recommended
% to change.
%
% thr2: threshold for extracting texture regions. Default value 12. Not
% recommended to change.
%
% Return Values
%
% band_vis_map: pixel-level banding visibility map.
%
% band_edge_list: cell array storing banding edges. Each cell is a banding edge.
%
% band_score: banding metric score. Larger means more banding.
%
% grad_mag: gradient magnitude feature map in case need further processing
%
% grad_dir: gradient direction feature map in case need further processing
%
% Required toolboxes: (should be included in path)
%
% drawedgelist.m
% edgelink.m
% filledgegaps.m
% findendsjunctions.m
% findisolatedpixels.m
% nonmaxsup.m
% circularstruct.m

%% check input and set default params
if ~exist('im_y','var') || isempty(im_y)
  error('Must have at least one non-empty input argument !');
end
if ~exist('imguidedfilt_wz','var') || isempty(imguidedfilt_ws)
  imguidedfilt_ws = 1;
end
if ~exist('thr1','var') || isempty(thr1)
  thr1 = 2;
end
if ~exist('thr2','var') || isempty(thr2)
  thr2 = 12;
end

% input image must be 1-channel
[~, ~, z] = size(im_y);
assert(z == 1, 'Input image must be 1-channel!');
% get image size
[im_height, im_width] = size(im_y);
% cvt double precesion
im_y = double(uint8(im_y));
% get #(pixels)
num_pixels = numel(im_y);
num_pixels_sqrt = sqrt(num_pixels);

%% Step1: Pre-processing
% Guided image filtering
if imguidedfilt_ws == 1
    im_y_filt = im_y;
else
    im_y_filt = imguidedfilter(im_y, ...
      'NeighborhoodSize', [imguidedfilt_ws imguidedfilt_ws]);
end

% MSCN feature extract in org image (MSCN Feature Map)
h = fspecial('gaussian', 9, 9/6);
im_mu = imfilter(im_y, h, 'replicate', 'same');
im_mu_sq = im_mu .* im_mu;
im_sigma = sqrt(abs(imfilter(im_y .* im_y, h, 'replicate', 'same') - im_mu_sq));
im_mscn = (im_y - im_mu) ./ (im_sigma + 1);

% Gradient feature map in *filtered* image using Sobel kernel (Grad Feature Map)
% initialize Sobel kernel
vx = [-1 0 +1];
g = [1 2 1];
fx = g' * vx;
fy = fx';
% convolution to get x- and y- gradient
grad_mag_x = imfilter(im_y_filt, fx, 'replicate', 'same');
grad_mag_y = imfilter(im_y_filt, fy, 'replicate', 'same');
% calculate gradient mag and dir
[grad_mag, grad_dir] = imgradient(grad_mag_x, grad_mag_y);

% calcualte Flat Pixel (FP), Text Pixel (TP), Band Pixel (BP) Feature Map

% double threshold followed by morphological majority filter
% Low threshold to get FP Feature Map
bmap_flat_pixels = grad_mag < thr1;
bmap_flat_pixels = medfilt2(bmap_flat_pixels, [3 3]);
% High threshold to get TP Feature Map
bmap_text_pixels = grad_mag > thr2;
bmap_text_pixels = medfilt2(bmap_text_pixels, [3 3]);
% Remainder is the BP Feature Map
bmap_band_pixels = ones(im_height, im_width) & ~bmap_text_pixels ...
                                             & ~bmap_flat_pixels;

%% Step2: banding edge detection

% set blk params
blksize_row = 9;
blksize_col = 9;
pad_row = floor(blksize_row/2);
pad_col = floor(blksize_col/2);

% pad image for border case
bmap_text_pixels_pad = padarray(bmap_text_pixels, [pad_row pad_col], ...
                            'replicate',  'both');

% Start from pre-defined BPs
[M, N] = find(bmap_band_pixels > 0);

% 1) Neighboring consistency check:
% BP points with any TP in neighbors are discarded
bmap_text_pixels_pad_integral = get_integral_image(bmap_text_pixels_pad);
for ii = 1:length(M)
  mm  = M(ii);
  nn  = N(ii);
  s_blk=integral_lookup(bmap_text_pixels_pad_integral,...
                        [mm, min(size(bmap_text_pixels_pad_integral, 1), mm+blksize_row -1)],...
                        [nn, min(size(bmap_text_pixels_pad_integral, 2), nn + blksize_col-1)]);
  if s_blk > 0
    bmap_band_pixels(mm, nn) = 0;
  end
end

% 2) Edge thining: Non-maxima suppression
grad_mag_band_pixels = grad_mag .* double(bmap_band_pixels);
% along the gradient normal orientation
k_radius = 1.5;
[grad_mag_nonmaxsup, ~] = nonmaxsup(grad_mag_band_pixels, ...
                                    abs(grad_dir), k_radius);

% refresh binary map
bmap_band_pixels = (bmap_band_pixels .* grad_mag_nonmaxsup) > 0;

% 3) Edge gap filling: fill small gaps in a binary edge map
k_fill_gap_size = 3;
bmap_band_pixels = filledgegaps(bmap_band_pixels, k_fill_gap_size);

% 4) Edge linking: linking pixels to form edges with 8-connectedneighbors.
try
  [band_edge_list, ~, ~] = edgelink(bmap_band_pixels);
catch
  % if didn't detect any edge, return none
  band_edge_list = {};
  band_vis_map = zeros(im_height, im_width);
  band_score = 0;
  return
end

%% Step3:  banding visibility estimation based on human vision model
% set blk params
blksize_row = 9;
blksize_col = 9;
pad_row = floor(blksize_row/2);
pad_col = floor(blksize_col/2);

% pad mscn image
im_mscn_pad = padarray(im_mscn, [pad_row pad_col], 'replicate', 'both');

% init banding visibility map
band_vis_map = zeros(im_height,im_width);

% set minimum edge length
k_min_edge_len = 10;

% traverse edge edge
for edge_idx = 1:length(band_edge_list)
  edge = band_edge_list{edge_idx};
  edge_len = size(edge,1);
  % discard edges shorter than k_min_edge_len pixels
  if edge_len < k_min_edge_len
    band_edge_list{edge_idx} = {};
    continue;
  end
  % edge cardinality weighting
  w_card_mask = vtf_card_mask(edge_len / num_pixels_sqrt);

  % traverse each pixel in current edge
  for pix_idx = 1:size(edge,1)
    mm = edge(pix_idx, 1);
    nn = edge(pix_idx, 2);

    blk_mscn = im_mscn_pad(mm:mm + blksize_row - 1, ...
                           nn:nn + blksize_col - 1);

    % banding edge contrast metric ( == gradient magnitude )
    band_edge_contrast = grad_mag(mm, nn);

    % Luminance masking
    local_lumi = im_mu(mm, nn);
    w_lumi_mask = vtf_lumi_mask( local_lumi );

    % Texture masking
    local_text = mean(mean(abs(blk_mscn)));
    w_text_mask = vtf_text_mask( local_text );

    % Overall banding score
    band_vis_map( mm, nn ) = ...
          band_edge_contrast * w_lumi_mask * w_text_mask * w_card_mask;
  end
end

%% Step4: Visual Importance Pooling
% average the most visible 80-percentile banding pixels
pct = 20;
banding_score_list = band_vis_map(band_vis_map > 0);
score_pct = prctile(banding_score_list, pct);
banding_score_list_pct = banding_score_list( banding_score_list > score_pct );
% get SI masking
w_si = vtf_si(get_si(im_y));
% BBAD_I score
band_score = w_si * mean(banding_score_list_pct);
end

% luminance masking visibility transfer function (VTF)
function weight = vtf_lumi_mask(l)
if l <= 0
  weight = 0;
elseif l <= 81
  % don't weighting down dark areas
  weight = 1;
elseif l <= 255
  % weighting down bright pixels
  weight = 1 - 1.6 * 10 ^ (-5) * (l - 81) ^ 2;
else
  weight = 0;
end
end
% texture masking VTF
function weight = vtf_text_mask(t)
if t <= 0.15
  weight = 1;
else
  weight = 1 / (( 1 + t - 0.15 ) ^ 5);
end
end
% edge cardinality masking VTF
function weight = vtf_card_mask(e)
weight = e ^ 0.5;
end
% get spatial information (SI) of Y
function si = get_si( Y )
[~, ~, z] = size(Y);
if z > 1
  Y = rgb2gray(Y);
end
% Cast to uint8 for compatibility with octave
[mag, ~] = imgradient(uint8(Y));
si = std2(mag);
end
% spatial information VTF
function weight = vtf_si(si)
weight = exp( -(si / 100) ^ 3 );
end
