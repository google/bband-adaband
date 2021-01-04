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
%
% ex, ey: the first(entry) pixel for the new segment
% map_segs: homosegment map
% map_segs_boundary: homosegment map for pixels on segments boundaries
% segs_entry: coodinates of the entry pixel for each segment
% im: 2D intensity image, visited pixels will be marked as -1
%
function [map_segs, map_segs_boundary, segs_entry] =...
    FindANewSeg(ex, ey, map_segs, map_segs_boundary, segs_entry, im)
  % The minimum segment size
  kMinSegNodes = 0.00001 * size(im, 1) * size(im, 2);

  % Four directions: up, left, bottom, right.
  NHOOD = [-1, 0; 0, -1; 1, 0; 0, 1];

  intensity = im(ex, ey);
  id = size(segs_entry, 1) + 1;

  num = 0;
  queue = [ex, ey];
  while size(queue, 1) ~= 0
    tx = queue(1, 1);
    ty = queue(1, 2);
    queue = queue(2:end, :);

    is_full_connected = true;
    for i = 1:size(NHOOD, 1)
      nx = tx + NHOOD(i, 1);
      ny = ty + NHOOD(i, 2);
      if (nx > 0 && nx <= size(im, 1) && ny > 0 && ny <= size(im, 2) ...
          && (im(nx, ny) == intensity ...
          || map_segs(nx, ny) == id ...
          || map_segs(nx, ny) == -id))
        if (im(nx, ny) == intensity && map_segs(nx, ny) == 0)
          map_segs(nx, ny) = -id;
          queue = [queue; nx, ny];
        end
      else
        is_full_connected = false;
      end
    end

%     im(tx, ty) = 0;
    if (is_full_connected == true)
      map_segs(tx, ty) = id;
      num = num + 1;
    end
  end

  if (num < kMinSegNodes)
    map_segs(map_segs == id | map_segs == -id) = 0;
  else
    [X, Y] = find(map_segs == id);
    for i = 1:length(X)
      tx = X(i);
      ty = Y(i);
      for j = 1:size(NHOOD, 1)
        nx = tx + NHOOD(j, 1);
        ny = ty + NHOOD(j, 2);
        if (nx > 0 && nx <= size(im, 1) && ny > 0 && ny <= size(im, 2) ...
            && ((im(nx, ny) == intensity && map_segs(nx, ny) == 0) ...
                || (map_segs(nx, ny) == -id)))
          map_segs_boundary(nx, ny) = id;
        end
      end
    end

    map_segs(map_segs == -id) = 0;

    [ex, ey] = find(map_segs_boundary == id, 1);
    segs_entry = [segs_entry; [ex, ey]];
  end
end
