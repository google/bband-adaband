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

function im_y_deband=video_postproc_deband_frame(im_y, im_y_pad_integral, band_vis_map, mask_proc_all, im_y_prev, im_y_deband_prev, map_half_ws, max_half_ws, min_edge_len, min_seg_area)

% get banding edge mask
band_vis_map_mask = band_vis_map > 0;

im_y_deband = im_y;

if ~isempty(im_y_prev),
   im_y_absframediff = abs( im_y - im_y_prev );
end

[vid_height, vid_width] = size(im_y);

%% 2) deband filter
% init dithering pattern just once for each frame
im_rand = 2 * (rand(vid_height,vid_width) - 0.5);
% filter all pixels in the mask_proc_all
[M_procpixs, N_procpixs] = find(mask_proc_all > 0);
for xxx = 1:length(M_procpixs)
  mmm = M_procpixs(xxx);
  nnn = N_procpixs(xxx);
  % randomize window size
  half_ws_tmp = randi(map_half_ws(mmm, nnn));
  if half_ws_tmp < 1, continue; end
  % extract local blk centered at current pixel
  M = [mmm + max_half_ws - half_ws_tmp, mmm + max_half_ws + half_ws_tmp];
  N = [nnn + max_half_ws - half_ws_tmp, nnn + max_half_ws + half_ws_tmp];
  blk_npixs = (2 * half_ws_tmp + 1) * (2 * half_ws_tmp + 1);
  blk_mu = integral_lookup(im_y_pad_integral, M, N) / blk_npixs;
  % dithering
  y_deband = blk_mu + im_rand(mmm, nnn);
  % apply spatial constrains
  if abs(y_deband - im_y(mmm, nnn)) <= 2
    if isempty(im_y_prev)
      im_y_deband(mmm, nnn) = y_deband;
    else
      % apply temporal constrains
      diffmag = abs(y_deband - im_y_deband_prev(mmm,nnn));
      diffsign = sign(y_deband - im_y_deband_prev(mmm,nnn));
      if diffmag <=  im_y_absframediff(mmm, nnn) + 1
        im_y_deband(mmm,nnn) = y_deband;
      else
        % soft thresholding larger change
        diffval = min(diffmag, im_y_absframediff(mmm, nnn ) + 1) * diffsign;
        im_y_deband(mmm,nnn) = im_y_deband_prev(mmm,nnn) + diffval;
      end
    end
  end
end

