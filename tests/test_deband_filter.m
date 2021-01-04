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

function test_deband_filter()
global save_expectations;
save_expectations = 0;

run_test_case('sky2_1s_trans_248_1920x1080', [1080, 1920], 1);
run_test_case('plant1_1s_trans_248_1920x1080', [1080, 1920], 1);
run_test_case('KristenAndSara_10frames_512x256', [256, 512], 10);

function run_test_case(name, size, num_frames)
global save_expectations;
mkdir('output')
rng(1011);
deband_filter_yuv420p(['testdata/', name, '.yuv420p'], ...
		      size, ['output/',name, '.yuv420p']);
failed = 0;
for i=1:num_frames,
  actual_y = read_yuv420_frame(['output/', name, '.yuv420p'], size, i-1);
  expected_y = read_yuv420_frame(['testdata/expectations/deband_filter/', ...
                                 name, '.yuv420p'], size, i-1);
  diff = actual_y - expected_y;
  if sum(abs(diff(:))) > 5,
    fprintf('%s has diff: %f\n', name, sum(abs(diff(:))));
    failed = 1;
  end
end
if failed
  fprintf('%s failed\n', name);
else
  fprintf('%s passed\n', name);
end
