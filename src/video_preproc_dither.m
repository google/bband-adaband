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

%% Scripts for applying pre-processing dithering to videos
%
%
% Instructions:
% Sctips to run proposed pre-processing dithering methods to input videos. Note
% that original videos are also passed to the rgb-yuv-rgb color conversion for
% fair comparison.
%
% change the parameters in %%Configs%% session and run this entire script.
%
% Required toolboxes: (should be included in path)
%
% 'ffmpeg'
% 'transcoder2' (google3)
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
curr_dir = pwd;
data_dir = 'data';
lib_dir = 'lib';
result_dir = 'result_preprocessing';
% video parameters
% if different framerates, should store in array [30,25,25,24..] and modify the
% call of video_framerate in Main Body to video_framerate(i)~
video_framerate = 30;
vid_width = 1280;
vid_height = 720;
% banding detection params
imguidedfilt_ws = 5;
% noise params
noise_level = 2;
half_ws = 16;
% transcoding setting: '248,247,137,136'...
output_formats = '247';
% output_formats = '248,247,137,136';
transcoder_abs_dir = '/google/src/cloud/zhengzhongtu/debanding_addNoise/google3/blaze-bin/video/transcoder2/tools/transcoder';
%%%%%%%%%%%%%%%

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

  %% pass original videos
  disp(['-- Processing ',name, '...']);
  % store extracted frames
  frames_raw_dir = fullfile(result_dir, 'frames_raw_dir');
  if ~exist(frames_raw_dir, 'dir')
    mkdir(frames_raw_dir);
  end
  % delete all png files if existed
  cd(frames_raw_dir)
  delete *.png
  cd(curr_dir)

  system(['ffmpeg -hide_banner -y -i ', video_name, ' ', ...
              fullfile(frames_raw_dir,'thumb%04d.png')]);
  num_frames = length(dir(frames_raw_dir))-2;

  frames_proc_dir = fullfile(result_dir, 'frames_proc_dir');
  if ~exist(frames_proc_dir, 'dir')
    mkdir(frames_proc_dir);
  end
  % delete all png files if existed
  cd(frames_proc_dir)
  delete *.png
  cd(curr_dir)

  % cvt each video for fair comparison
  for jj = 1:num_frames
    im_rgb = im2double(imread(fullfile(frames_raw_dir, ...
      sprintf('thumb%04d.png',jj))));
    % cvt to yuv
    im_yuv = rgb2ycbcr(im_rgb);
    % only process y plane
    im_rgb = ycbcr2rgb(im_yuv);
    imwrite(im2uint8(im_rgb), fullfile(frames_proc_dir, ...
        [name, sprintf('_frame%04d.png',jj)]));
  end
  % video output name
  out_dir = fullfile(result_dir, name);
  if ~exist(out_dir, 'dir')
    mkdir(out_dir);
  end
  video_raw_name = fullfile(out_dir, [name, '.y4m']);
  % concate frames to video
  cmd = ['ffmpeg -hide_banner -y', ' -framerate ', ...
    num2str(video_framerate),' -i ', fullfile(frames_proc_dir, ...
    [name, '_frame%04d.png']), ' -pix_fmt yuv420p', ' -start_number 1', ...
    ' -s ', num2str(vid_width), 'x',num2str(vid_height), ' ', video_raw_name];
  system(cmd);

  %% Banding detection and adding noise
  % extracted frames from videos
  frames_raw_dir = fullfile(result_dir, 'frames_raw_dir');
  if ~exist(frames_raw_dir, 'dir')
    mkdir(frames_raw_dir);
  end
  % delete all png files if existed
  cd(frames_raw_dir)
  delete *.png
  cd(curr_dir)

  system(['ffmpeg -hide_banner -y -i ', video_name, ' ', ...
              fullfile(frames_raw_dir,'thumb%04d.png')]);
  num_frames = length(dir(frames_raw_dir)) - 2;

  frames_proc_dir = fullfile(result_dir, 'frames_proc_dir');
  if ~exist(frames_proc_dir, 'dir')
    mkdir(frames_proc_dir);
  end
  % delete all png files if existed
  cd(frames_proc_dir)
  delete *.png
  cd(curr_dir)

  banding_scores_BBAD_tmp  = zeros(1, num_frames);
  banding_score_BBAD_I_frames{1,i} = zeros(1, num_frames);
  % store previous frame
  im_prev_frame = [];
  k = 0;

  for jj = 1:num_frames
    k = k + 1;
    im_rgb = im2double(imread(fullfile(frames_raw_dir, ...
                          sprintf('thumb%04d.png',jj))));
    % cvt to yuv in double precision
    im_yuv = rgb2ycbcr(im_rgb);
    % only process y plane
    im_yyy = im_yuv(:,:,1);

    % Banding region detectino by BBAD
    disp('-- Running BBAD_I')
    tic
    [band_vis_map, band_edge_list, band_score, grad_mag, grad_dir] = ...
                        BBAD_I(im2uint8(im_yyy), imguidedfilt_ws);
    toc
    %% Adding noise
    frames_proc_dir = fullfile(result_dir, 'frames_proc_dir');
    if ~exist(frames_proc_dir, 'dir')
      mkdir(frames_proc_dir);
    end
    % delete all png files if existed

    disp('-- Pre-processing by dithering ...')
    % Gray scale processing -- only processing Y-ch
    % initialize local processsing
    im_yyy_local_dith = im_yyy;
    tic
    [M_bandpixs,N_bandpixs] = find(band_vis_map > 0);
    % only process non-textured pixels
    mask_proc = grad_mag <= 12;

    % apply uniform noise around each pixel
    for ii = 1:length(M_bandpixs)
      mm = M_bandpixs(ii);
      nn = N_bandpixs(ii);
      blk_mm_range = max([1, mm - half_ws]): ...
                     min([mm + half_ws, vid_height]);
      blk_nn_range = max([1, nn - half_ws]): ...
                     min([nn + half_ws, vid_width ]);
      im_yyy_local_dith(blk_mm_range, blk_nn_range) = ...
          im_yyy( blk_mm_range, blk_nn_range ) + ...
          2*noise_level.*(rand(length(blk_mm_range), ...
                                length(blk_nn_range))-0.5)./255;
    end
    % replace
    im_yyy( mask_proc ) = im_yyy_local_dith( mask_proc );

    toc
    % writing our dithered frames.
    im_yuv(:,:,1) = im_yyy;
    im_rgb = ycbcr2rgb(im_yuv);
    imwrite(im2uint8(im_rgb), fullfile(frames_proc_dir,[name, '_dither', ...
      sprintf('_frame%04d.png',jj)]));

    %% modulated by temporal score
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

  %% concate frames to videos

  out_dir = fullfile(result_dir,'pre_proc_videos', name);
  if ~exist(out_dir,'dir')
    mkdir(out_dir);
  end
  video_name_out = fullfile(out_dir, [name, '_dither', '.y4m']);
  cmd = ['ffmpeg -hide_banner -y', ' -framerate ', num2str(video_framerate),...
    ' -i ', fullfile(frames_proc_dir,[name,'_dither','_frame%04d.png']), ...
      ' -pix_fmt yuv420p', ' -start_number 1', ' -s ', ...
      num2str(vid_width), 'x',num2str(vid_height), ' ', video_name_out];
  system(cmd);
  % remove the tmp frames
  cd(frames_proc_dir)
  delete *.png
  cd(curr_dir)

  %% trancode videos to 248,247...
  out_dir = fullfile(result_dir, 'pre_proc_videos', name);
  if ~exist(out_dir,'dir'),  mkdir(out_dir); end
  video_name_out = fullfile(out_dir, [name, '_dither', '.y4m']);
  video_in_abs_path = fullfile(curr_dir, video_name_out);
  video_trans_abs_path = fullfile(out_dir, [name, '_dither']);

  cmd = [transcoder_abs_dir, ...
  ' --input_file=', video_in_abs_path, ' --output_file_prefix=', ...
  video_trans_abs_path, ' --output_formats=', output_formats, ...
  ' --logtostderr'];
  system(cmd);
end



