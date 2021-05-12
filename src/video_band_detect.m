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

%% Scripts for evaluating banding detection benchmarks on videos
%
% Instructions:
% Sctips to run proposed video level as well as frame level banding detector on
% given videos stored in a directory. It will automaticly run detection on all
% the videos stored within this dir.
%
% change the parameters in %%Configs%% session and run this entire script.
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
% FCDR.m
% BBAD_I.m

close all;
clear;
warning('off','all');

%% Configs %%
curr_dir = pwd;
data_dir = 'data';
result_dir = 'result_banding_detect';
lib_dir = 'lib';
% video parameters
% set frame sampling step (15 means 1 out of every 15 frames)
sampling_step = 15;
% choose which metric to evaluate: 'BBAD', 'FCDR'
eval_flags = { 'BBAD' };
% eval_flags = { 'BBAD', 'FCDR' };
%%%%%%%%%%%%%

%% Main body
% init
if ~exist(result_dir, 'dir')
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
json_BBAD_V = cell(1, length(test_seqs));

% frame level banding score
banding_score_BBAD_I_frames = cell(1, length(test_seqs));

for i = 1:length(test_seqs)
    video_name = fullfile(data_dir, test_seqs{i});
    [filepath, name, ext] = fileparts(video_name);
    disp(['-- Processing ', video_name, '...']);

    % extracted frames from videos
    frames_raw_dir = fullfile(result_dir, 'frames_raw_dir');
    if ~exist(frames_raw_dir,'dir')
      mkdir(frames_raw_dir);
    end
    % delete all png files if existed
    cd(frames_raw_dir)
    delete *.png
    cd(curr_dir)
    system(['ffmpeg -hide_banner -y -i ', video_name, ' ', ...
                fullfile(frames_raw_dir,'thumb%04d.png')]);
    num_frames = length(dir(frames_raw_dir)) - 2;
    % banding detect data
    banding_score_BBAD_I_frames{1,i} = [];
    im_prev_frame = [];
    k = 0;
    % read one frame
    for jj = 1:num_frames
        k = k+1;
        im = imread(fullfile(frames_raw_dir, ...
                              sprintf('thumb%04d.png',jj)));
        % frame sampling
        if mod(k,sampling_step)~=0
            im_prev_frame = im;
            continue;
        end
        % color conversion
        im_yuv = double( im2uint8(rgb2ycbcr(im2double(im))));
        im_y = im_yuv(:,:,1);
        % evaluate FCDR
        if find(ismember(eval_flags, 'FCDR'))
            algo = 'FCDR';
            disp('-- Running FCDR ')
            tic;
            [fcc_map_1, fcc_map_2, fcc_map_3] = FCDR(im);
            toc;
            % write results
            out_dir = fullfile(result_dir, 'FCDR_output', name);
            if ~exist(out_dir, 'dir')
              mkdir(out_dir);
            end
            mat_name_out = fullfile(out_dir, sprintf('frame%04d_FCDR.mat',k));
            save(mat_name_out, 'fcc_map_1', 'fcc_map_2', 'fcc_map_3');
        end

        % Evaluate BBAD
        if find(ismember(eval_flags, 'BBAD'))
            algo = 'BBAD';
            disp('-- Running BBAD')
            tic;
            [band_vis_map, band_edge_list, band_score, grad_mag, grad_dir] = ...
                            BBAD_I(im_y, 1);
            toc;
            % modulated by temporal score
            if isempty(im_prev_frame)
                ti = 0;
            else
                ti = std2(rgb2gray(uint8(im)) - rgb2gray(uint8(im_prev_frame)));
            end
            % get TI masking weight
            w_ti = exp( - (ti / 20) ^ 2 ) ;
            % calculate BBAD_I score for curr frame
            banding_score_BBAD_I_frames{1,i}(end+1) = w_ti * band_score;
            % write results
            out_dir = fullfile(result_dir, 'BBAD_output', name);
            if ~exist(out_dir, 'dir')
              mkdir(out_dir);
            end
            mat_name_out = fullfile(out_dir, sprintf('frame%04d_BBAD.mat',k));
            save(mat_name_out, 'band_vis_map', 'band_edge_list', ...
                 'band_score', 'grad_mag', 'grad_dir');
        end
        im_prev_frame = im;
    end
    banding_score_BBAD_V(1,i) = mean(banding_score_BBAD_I_frames{1,i});
    json_BBAD_V{1,i} = struct('input_file', name, 'BBAND', banding_score_BBAD_V(1,i));
end
s.all_scores = json_BBAD_V;
fileID = fopen('bband_results.json', 'w');
fprintf(fileID, jsonencode(s));
fclose(fileID);

