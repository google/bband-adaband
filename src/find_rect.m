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

function [I,J] = find_rect(map_segs, min_x, max_x, a, b, n)
% Finds pixels in map_segs that are equal to either a or b looking
% only within the rectangle defined by min_x to max_x.
rect = map_segs(min_x(1):max_x(1), min_x(2):max_x(2));
if nargin < 6
   [I, J] = find((rect == a) | (rect == b));
else
   [I, J] = find((rect == a) | (rect == b), n);
end
I = I + min_x(1) - 1;
J = J + min_x(2) - 1;