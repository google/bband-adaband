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

function r=integral_lookup(im, i, j)
% Perform a look up on the integral image of the region defined by
% i(1):i(2) and j(1):j(2).
% im should be an integral image.

top = 0;
left = 0;
top_left = 0;
if i(1) > 1,
  top = im(i(1) - 1, j(2));
end
if j(1) > 1,
  left = im(i(2), j(1) - 1);
end
if (i(1) > 1) & (j(1) > 1)
  top_left = im(i(1) - 1, j(1) - 1);
end
r=im(i(2), j(2)) - top - left + top_left;
