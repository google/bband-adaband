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

%% Scripts for applying post-processing debanding filter to videos
%
%
% Instructions:
% Sctips to run proposed post-processing debanding filter to input videos. This
% filtering is only applied to the Y-channel.
%
% change the parameters in %%Configs%% session and run this entire script.
%
% Required toolboxes: (should be included in path)
%
% 'ffmpeg'
% drawedgelist.m
% edgelink.m
% filledgegaps.m
% findendsjunctions.m
% findisolatedpixels.m
% nonmaxsup.m
% circularstruct.m
% BBAD_I.m

close all;
clear;
warning('off','all');

%% Configs %%
ffmpeg = 'ffmpeg'
curr_dir = pwd;
data_dir = 'data';
lib_dir = 'lib';
result_dir = 'result_postprocessing';
% video parameters
% if different framerates, should store in array [30,25,25,24..] and modify the
% call of video_framerate in Main Body to video_framerate(i)~
vid_framerate = 30;
% banding detection params
imguidedfilt_ws = 1;
% deband filter params
% maximum half window size
max_half_ws = 64;
% edge length threshold, skip edges shorter than min_edge_len
min_edge_len = 16;
% segment area threshold, skip segs smaller than this valu
min_seg_area = 32*32;
%%%%%%%%%%%%%

%% Main Body
if ~exist(result_dir,'dir')
  mkdir(result_dir);
end
addpath(genpath(lib_dir))

% parse data names
test_seqs = dir(data_dir);
test_seqs = {test_seqs.name};
% dir results always like  '.','..','video1.mp4','video2.mp4',...
% so we start from 3 to exclude '.','..'
test_seqs = test_seqs(3:end);

% video level banding score
banding_score_BBAD_V = zeros(1, length(test_seqs));
% frame level banding score
banding_score_BBAD_I_frames = cell(1, length(test_seqs));

% for each sequence
for i = 1:length(test_seqs)
  video_name = fullfile(data_dir, test_seqs{i});
  [filepath, name, ext] = fileparts(video_name);

  disp(['-- Processing ',name, '...']);

  %% Banding detection and debanding filter
  % Here we can do banding detection and deband filtering together, but
  % instead, we do it separately for easier debug.
  % extracted frames from videos
  frames_raw_dir = fullfile(result_dir, 'frames_raw_dir');
  if ~exist(frames_raw_dir,'dir')
    mkdir(frames_raw_dir);
  end
  % delete all png files if existed
  cd(frames_raw_dir)
  delete *.png
  cd(curr_dir)

  system([ffmpeg, ' -hide_banner -y -i ', video_name, ' ', ...
              fullfile(frames_raw_dir,'thumb%04d.png')]);
  num_frames = length(dir(frames_raw_dir)) - 2;
  % banding detection data
  banding_score_BBAD_I_frames{1,i} = zeros(1, num_frames);
  % store previous frame
  im_prev_frame = [];
  k = 0;

  for jj = 1:num_frames
      k = k + 1;

      im_rgb = imread(fullfile(frames_raw_dir, ...
                            sprintf('thumb%04d.png',jj)));
      % cvt to yuv in double precision
      im_yuv = double(rgb2ycbcr(im_rgb));
      % only process y plane
      im_yyy = im_yuv(:,:,1);

      % Potential Banding region detectino by BBAD
      disp(sprintf('-- Running BBAD_I %d', jj))
      tic
      [band_vis_map, band_edge_list, band_score, grad_mag, grad_dir] = ...
                          BBAD_I(im_yyy, imguidedfilt_ws);
      toc
      % modulated by temporal score
      if isempty(im_prev_frame)
          ti = 0;
      else
          ti = std2(rgb2gray(uint8(im_rgb)) - rgb2gray(uint8(im_prev_frame)));
      end
      % get TI masking weight
      w_ti = exp( - (ti / 20) ^ 2 ) ;
      % calculate BBAD_I score for curr frame
      banding_score_BBAD_I_frames{1,i}(jj) = w_ti * band_score;
      % write results
      out_dir = fullfile(result_dir, 'BBAD_output', name);
      if ~exist(out_dir,'dir')
        mkdir(out_dir);
      end
      mat_name_out = fullfile(out_dir, sprintf('frame%04d_BBAD.mat',k));
      save(mat_name_out, 'band_vis_map', 'band_edge_list', ...
               'band_score', 'grad_mag', 'grad_dir');
      im_prev_frame = im_rgb;
  end
  banding_score_BBAD_V(1,i) = mean(banding_score_BBAD_I_frames{1,i});

  %% Proposed adaptive average debanding filter
  % notice: this module must be executed after previous module finished
  % extract frames
  frames_raw_dir = fullfile(result_dir, 'frames_raw_dir');
  if ~exist(frames_raw_dir, 'dir')
    mkdir(frames_raw_dir);
  end
  % delete all png files if existed
  cd(frames_raw_dir)
  delete *.png
  cd(curr_dir)
  system([ffmpeg, ' -hide_banner -y -i ', video_name, ' ', ...
              fullfile(frames_raw_dir, 'thumb%04d.png')]);

%     frames_deband_dir = fullfile(result_dir, 'frames_deband_dir');
%     if ~exist(frames_deband_dir, 'dir')
%       mkdir(frames_deband_dir); 
%     end
  k_method = 'proposed_adap_deband_filt';
  num_frames = length(dir(frames_raw_dir)) -2;
  % for each frame apply deband filter
  for jj = 1:num_frames
      fprintf('--Processing %04d th frame...\n', jj);
      % Load pre-processed banding visibility map from previous module
      disp('-- Loading BBAD results...')
      % if didn't load, skip this frame
      try
          out_dir = fullfile(result_dir, 'BBAD_output', name);
          mat_name_out = fullfile(out_dir, sprintf('frame%04d_BBAD.mat',jj));
          load(mat_name_out);
      catch
          continue;
      end

      if jj > 1
        % get prev orig frame
        im_rgb_prev = imread( fullfile(frames_raw_dir, ...
                                sprintf('thumb%04d.png',jj-1)));
        im_yuv_prev = double(im2uint8(rgb2ycbcr(im2double(im_rgb))));
        im_y_prev   = im_yuv(:,:,1);
        % get prev deband frame
        deband_result_dir = fullfile(result_dir, 'video_deband_results', name);
        if ~exist(deband_result_dir,'dir')
          mkdir(deband_result_dir);
        end
        im_rgb_deband_prev = imread(fullfile(deband_result_dir, ...
                  [name, '_', k_method, sprintf('_frame%04d.png', jj-1)]));
        im_yuv_deband_prev = double(im2uint8(rgb2ycbcr(im2double(im_rgb_deband_prev))));
        im_y_deband_prev = im_yuv_deband_prev(:,:,1);
      else
        im_y_prev = [];
        im_y_deband_prev = [];
      end
      % get current frame
      im_rgb = imread(fullfile(frames_raw_dir, sprintf('thumb%04d.png', jj)));
      % cvt to yuv
      im_yuv = double(im2uint8(rgb2ycbcr(im2double(im_rgb))));
      % only process y plane
      im_y = im_yuv(:,:,1);

      [mask_proc_all, im_y_pad_integral, map_half_ws] = video_postproc_deband_get_mask(...
          im_y, grad_mag, grad_dir, band_edge_list, max_half_ws, min_edge_len, min_seg_area);
      im_y_deband=video_postproc_deband_frame(im_y, im_y_pad_integral, band_vis_map, mask_proc_all, ...
	 im_y_prev, im_y_deband_prev, map_half_ws, max_half_ws, min_edge_len, min_seg_area);

    im_yuv(:,:,1) = im_y_deband;
    im_rgb_deband = ycbcr2rgb(uint8(im_yuv));

    deband_result_dir = fullfile(result_dir, 'video_deband_results', name);
    if ~exist(deband_result_dir,'dir')
      mkdir(deband_result_dir);
    end
    im_name_out = fullfile(deband_result_dir, [name,'_',k_method, sprintf('_frame%04d.png', jj)]);
    imwrite(im_rgb_deband, im_name_out);
  end

  %% concate frames to videos
  deband_result_dir = fullfile(result_dir,'video_deband_results', name);
  if ~exist(deband_result_dir,'dir'),  mkdir(deband_result_dir); end
  video_name_out = fullfile(deband_result_dir, [name, '_', k_method '.y4m']);
  [video_height, video_width, ~] = size(im_yuv);
  cmd = [ffmpeg, ' -hide_banner -y', ' -framerate ', ...
    num2str(round(vid_framerate)), ' -i ', ...
    fullfile(deband_result_dir ,[name,'_',k_method,'_frame%04d.png']), ...
      ' -pix_fmt yuv420p', ' -start_number 1', ...
      ' -s ', num2str(vid_width), 'x',num2str(vid_height), ' ', video_name_out, '>/dev/null'];
  system(cmd);
end
