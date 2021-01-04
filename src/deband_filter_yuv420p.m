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

function deband_filter_yuv420(input_file, input_size, output_file)
%% Apply debanding filter to YUV420P input and write to output..
%
% Input Arguments
%
% input_file: Path to input yuv420p raw video.
% input_size: Dimensions of input frame in matlab format [height, width]
% output_file: Path to output yuv420p raw video.

imguidedfilt_ws = 1;
max_half_ws = 64;
min_edge_len = 16;
min_seg_area = 32*32;

im_y_prev = [];
im_y_deband_prev = [];
output_fid = fopen(output_file, "w");

num_frames=1000;
for f=1:num_frames,
  [im_y, u, v] = read_yuv420_frame(input_file, input_size, f - 1);
  if isempty(im_y)
    break
  end
  [band_vis_map, band_edge_list, band_score, grad_mag, grad_dir] = BBAD_I(im_y, imguidedfilt_ws);
  [mask_proc_all, im_y_pad_integral, map_half_ws] = video_postproc_deband_get_mask(...
    im_y, grad_mag, grad_dir, band_edge_list, max_half_ws, min_edge_len, min_seg_area);

  im_y_deband=video_postproc_deband_frame(im_y, im_y_pad_integral, band_vis_map, mask_proc_all, ...
    im_y_prev, im_y_deband_prev, map_half_ws, max_half_ws, min_edge_len, min_seg_area);

  yuv420 = {im_y_deband, u, v};
  FileUtils.write420(output_fid, yuv420);

  im_y_prev = im_y;
  im_y_deband_prev = im_y_deband;
end

fclose(output_fid);

