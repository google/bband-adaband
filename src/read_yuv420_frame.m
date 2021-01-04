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

function [y, u, v] = read_yuv420_frame(input_yuv, input_size, frame_number)
% Reads a YUV420P frame from the input video.
%
% If frame_number is beyond the file end, return empty arrays.
if nargin < 3
   frame_number = 0;
end
f = fopen(input_yuv, 'rb');
frame_size = input_size(1) * input_size(2) * 3 / 2;
if fseek(f, frame_number * frame_size, -1) < 0,
  y = [];
  u = [];
  v = [];
  fclose(f);
  return
end
y = fread(f, input_size(2:-1:1), 'uint8')';
u = fread(f, input_size(2:-1:1) / 2, 'uint8')';
v= fread(f, input_size(2:-1:1) / 2, 'uint8')';
fclose(f);