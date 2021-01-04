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

function [mask_proc_all, im_y_pad_integral, map_half_ws]=video_postproc_deband_get_mask(im_y, grad_mag, grad_dir, band_edge_list, max_half_ws, min_edge_len, min_seg_area)

% The following globals are reused between calls to FindANewSeg so that
% we don't have to allocate big images thousands of time per call
global map_segs;
global map_segs_boundary;
global seq_queue;
im_y_pad = padarray( im_y, [max_half_ws max_half_ws], ...
                     'replicate', 'both');
[im_y_pad_integral, im_y_pad_integral2] = get_integral_image(im_y_pad);

grad_mag_pad = uint16(padarray(grad_mag, [max_half_ws max_half_ws], ...
                               'replicate','both'));

seg_queue = zeros(size(im_y, 1) * size(im_y, 2), 2);
map_segs = zeros(size(im_y));
map_segs_boundary = zeros(size(im_y));

% mask to record all the processed pixels
%mask_proc_all = false(size(im_y));
% half window size map. updated for each loop
map_half_ws = zeros(size(im_y));
% connected component map, if a segment is already searched, we labeled
% it to prevent searching again.
map_conn_comp = zeros( size(im_y) );
% init segment label
idx_conn_comp = 1;
conn_comp = cell(length(band_edge_list), 1);
ran = zeros(length(band_edge_list), 1);
%% 1) get adaptive half_ws for each banding pixel
% process every edge independently
for ii = 1:length(band_edge_list)
  % for debug, print the %process
  if mod(ii, floor(0.1*length(band_edge_list)))==0
    disp([num2str(ii), ', ', ...
          num2str(ii/length(band_edge_list)*100), '%']);
  end
  % get current edge
  edge = band_edge_list{ii};
  edge_length = size(edge,1);
  % ignore short edges
  if edge_length < min_edge_len
    continue;
  end
  % using the midpoint of edge as seed to find conn component
  mid_point = edge(floor(edge_length/2), :);
  % gradient direction
  grad_angle = grad_dir(mid_point(1), mid_point(2));
  % 8-conn neighbors search
  xy_offset = (grad_angle >= -22.5 && grad_angle < 22.5) * [0 1] + ...
              (grad_angle >= 22.5 && grad_angle < 67.5) * [-1 1] + ...
              (grad_angle >= 67.5 && grad_angle < 112.5) * [-1 0] + ...
              (grad_angle >= 112.5 && grad_angle < 157.5) * [-1 -1] + ...
              (grad_angle >= 157.5 && grad_angle < -157.5) * [0 -1] + ...
              (grad_angle >= -157.5 && grad_angle <- 112.5) * [1 -1] + ...
              (grad_angle >= -112.5 && grad_angle < -67.5) * [1 0] + ...
              (grad_angle >= -67.5 && grad_angle < -22.5) * [1 1];
  % get seeds for episeg and hypseg around the edge
  entry_episeg = mid_point + xy_offset;
  % if current seg hasn't been visited, do connected component label
  if map_conn_comp(entry_episeg(1), entry_episeg(2)) == 0
    % 8-conn comp labeling to get the segs
    [min_x, max_x] = FindANewSeg2(entry_episeg(1), entry_episeg(2), idx_conn_comp, im_y);
    [M_procpixs_epi, N_procpixs_epi] = find_rect(map_segs, min_x, max_x, idx_conn_comp, idx_conn_comp);
    map_segs(sub2ind(size(map_segs), M_procpixs_epi, N_procpixs_epi)) = 0;
    
    map_conn_comp(sub2ind(size(map_conn_comp), M_procpixs_epi, N_procpixs_epi)) = idx_conn_comp;
    conn_comp{idx_conn_comp} = [M_procpixs_epi, N_procpixs_epi];
    idx_conn_comp = idx_conn_comp + 1;
  else
    M_procpixs_epi = conn_comp{map_conn_comp(entry_episeg(1),entry_episeg(2))}(:, 1);
    N_procpixs_epi = conn_comp{map_conn_comp(entry_episeg(1),entry_episeg(2))}(:, 2);
  end
  % get width and append
  area_episeg = length(M_procpixs_epi);
  if area_episeg < min_seg_area
    % if area is small, it's likely a block, so assign
    % sqrt(area) to it,(optimized for potential block segments)
    d_episeg = sqrt(area_episeg);
  else
    d_episeg = min(2 * max_half_ws, area_episeg/edge_length);
  end

  % do the same for hypseg
  entry_hypseg = mid_point - xy_offset;
  if map_conn_comp(entry_hypseg(1), entry_hypseg(2)) == 0
    [min_x, max_x] = FindANewSeg2(entry_hypseg(1), entry_hypseg(2), idx_conn_comp, im_y);
    [M_procpixs_hyp, N_procpixs_hyp] = find_rect(map_segs, min_x, max_x,idx_conn_comp, idx_conn_comp);
    map_segs(sub2ind(size(map_segs), M_procpixs_hyp, N_procpixs_hyp)) = 0;
    map_conn_comp(sub2ind(size(map_conn_comp), M_procpixs_hyp, N_procpixs_hyp)) = idx_conn_comp;
    conn_comp{idx_conn_comp} = [M_procpixs_hyp, N_procpixs_hyp];
    idx_conn_comp = idx_conn_comp + 1;
  else
    M_procpixs_hyp = conn_comp{map_conn_comp(entry_hypseg(1),entry_hypseg(2))}(:, 1);
    N_procpixs_hyp = conn_comp{map_conn_comp(entry_hypseg(1),entry_hypseg(2))}(:, 2);
  end
  area_hypseg = length(M_procpixs_hyp);
  if area_hypseg < min_seg_area
    d_hypseg = sqrt(area_hypseg);
  else
    d_hypseg = min(2*max_half_ws, area_hypseg/edge_length);
  end

  % throw very narrow seg's window size
  if d_episeg < 4 && d_hypseg < 4, continue; end
  % throw very unbalanced seg's window size
  if abs(d_episeg-d_hypseg)/(d_episeg+d_hypseg) > 0.95, continue; end

  % proc pixs in episeg
  epi_conn = map_conn_comp(entry_episeg(1),entry_episeg(2));
  if (epi_conn == 0) | (ran(epi_conn) ~= d_episeg),
    half_ws_episeg = find_adaptive_window_size(im_y, im_y_pad, im_y_pad_integral, ...
					       im_y_pad_integral2, d_episeg, ...
					       M_procpixs_epi, N_procpixs_epi, max_half_ws);
    I = sub2ind(size(map_half_ws), M_procpixs_epi, N_procpixs_epi);
    map_half_ws(I) = max(half_ws_episeg, map_half_ws(I));
    if epi_conn > 0,
      ran(epi_conn) = d_episeg;
    end
  end

  % do the same for hypseg
  hyp_conn = map_conn_comp(entry_hypseg(1),entry_hypseg(2));
  if (hyp_conn == 0) | (ran(hyp_conn) ~= d_hypseg),
    half_ws_hypseg = find_adaptive_window_size(im_y, im_y_pad, im_y_pad_integral, ...
					       im_y_pad_integral2, d_hypseg, ...
					       M_procpixs_hyp, N_procpixs_hyp, max_half_ws);
    I  = sub2ind(size(map_half_ws), M_procpixs_hyp, N_procpixs_hyp);
    map_half_ws(I) = max(half_ws_hypseg, map_half_ws(I));
    if (hyp_conn > 0)
      ran(hyp_conn) = d_hypseg;
    end
  end
end
% update mask of pixels processed (assume each processed pixel has map_half_ws set)
mask_proc_all = map_half_ws > 0;

% Finds a new homo segment
%
% ex, ey: the first(entry) pixel for the new segment
% map_segs: homosegment map
% map_segs_boundary: homosegment map for pixels on segments boundaries
% segs_entry: coodinates of the entry pixel for each segment
% im: 2D intensity image, visited pixels will be marked as -1
%
function [min_x, max_x] = FindANewSeg2(ex, ey, id, im)
  global seg_queue;
  global map_segs;
  global map_segs_boundary;
  % The minimum segment size
  kMinSegNodes = 0.00001 * size(im, 1) * size(im, 2);

  % Four directions: up, left, bottom, right.
  NHOOD = [-1, 0; 0, -1; 1, 0; 0, 1];

  intensity = im(ex, ey);

  num = 0;
  nq = 1;
  seg_queue(nq, :) = [ex, ey];
  min_x = [ex, ey];
  max_x = [ex, ey];
  while nq > 0,
    tx = seg_queue(nq, 1);
    ty = seg_queue(nq, 2);
    nq = nq - 1;

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
	  nq = nq + 1;
          seg_queue(nq, :) = [nx, ny];
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
    min_x = min(min_x, [tx, ty]);
    max_x = max(max_x, [tx, ty]);
  end

  min_x = max(min_x - 2, [1, 1]);
  max_x = min(max_x + 2, [size(im, 1), size(im, 2)]);
  if (num < kMinSegNodes)
    [I,J] = find_rect(map_segs, min_x, max_x, id, -id);
    map_segs(sub2ind(size(map_segs), I, J)) = 0;
  else
    [X, Y] = find_rect(map_segs, min_x, max_x, id, id);
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

    [I,J] = find_rect(map_segs, min_x, max_x, -id, -id);
    map_segs(sub2ind(size(map_segs), I, J)) = 0;

    [ex, ey] = find_rect(map_segs_boundary, min_x, max_x, id, id);
    % Union the boundary into the map segs
    map_segs(sub2ind(size(map_segs), ex, ey)) = id;

    % Clear the boundary too so it can be re-used for next call
    map_segs_boundary(sub2ind(size(map_segs_boundary), ex, ey)) = 0;
  end
