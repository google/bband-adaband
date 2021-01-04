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

function test_BBAD_I()
global save_expectations;
save_expectations = 0;

run_test_case('sky2_1s_trans_248_1920x1080', [1080, 1920])
run_test_case('plant1_1s_trans_248_1920x1080', [1080, 1920])

function run_test_case(name, size)
global save_expectations;

[y,u,v] = read_yuv420_frame(['testdata/', name, '.yuv420p'], size);
start_time=tic;
[actual_band_vis_map, actual_band_edge_list, actual_band_score, ~, ~] = BBAD_I(y, 1);
subplot(1, 2, 1);
imshow(y);
subplot(1, 2, 2);
imshow(actual_band_vis_map);
if save_expectations
  band_vis_map = actual_band_vis_map;
  band_score = actual_band_score;
  save(['testdata/expectations/BBAD_I/', name, '.mat'], ...
       'band_vis_map', 'band_score');
else
  load(['testdata/expectations/BBAD_I/', name, '.mat']);
  if sum(abs(actual_band_vis_map(:) - band_vis_map(:))) > 1
    fprintf('%s: Band vis map is different\n', name)
  else
    fprintf('%s: pass\n', name)
  end
end
toc(start_time)

