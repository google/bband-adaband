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

function window_size = find_adaptive_window_size(im_y, im_y_pad, im_y_pad_integral, im_y_pad_integral2, d, M, N, max_half_ws)
% Search for the adaptive window size by looking for the biggest window that
% has variance below a threshold.
%
% Input arguments:
% im_y: input image
% im_y_padded: padded input image (by half_ws_adapted)
% M, N: Pixel indices to do the search on.
%
% Output arguments:
% Outputs a map that has the window size for each pixel in [M, N]

window_size = zeros(length(M), 1);
for jjj = 1:length(M)
  mmm = M(jjj);
  nnn = N(jjj);
  half_ws_adapted = floor(d / 2);
  % make it adaptive: banding edge's local window should be low
  % variance. If local variance exceeds 4, half_ws scaled down by
  % 2 until local variance are below 4.
  while(half_ws_adapted > 1)
    blk_mmm_range = [mmm + max_half_ws - half_ws_adapted, mmm + max_half_ws + half_ws_adapted];
    blk_nnn_range = [nnn + max_half_ws - half_ws_adapted, nnn + max_half_ws + half_ws_adapted];
    y_sum = integral_lookup(im_y_pad_integral, blk_mmm_range, blk_nnn_range);
    y_sum2 = integral_lookup(im_y_pad_integral2, blk_mmm_range, blk_nnn_range);
    npix = (blk_mmm_range(2) - blk_mmm_range(1) + 1) * (blk_nnn_range(2) - blk_nnn_range(1) + 1);

    % Used unbiased estimate of variance (E(X^2)-E(X)^2)*npix/(npix-1)
    v = (y_sum2 / (npix - 1) - (y_sum / npix) * (y_sum / (npix - 1)));

    % Any differences between original implementation and the integral image is likely
    % due to slight differences in variance computation. Can debug here by enabling
    % the following check.
    if 0,
      y_blk = im_y_pad( blk_mmm_range(1):blk_mmm_range(2), blk_nnn_range(1):blk_nnn_range(2) );
      v2 = var(y_blk(:));
      if (v2 <= 4) ~= (v <= 4),
        fprintf("%f %f %f %f\n", v, v2, y_sum/npix, mean(y_blk(:)))
      end
    end

    if v <= 4 || half_ws_adapted <= 1, break;
    else
      half_ws_adapted = floor(half_ws_adapted / 2);
    end
  end
  window_size(jjj) = half_ws_adapted;
end
