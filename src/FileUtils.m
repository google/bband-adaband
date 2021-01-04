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

% Force matlab-mode in emacs: -*-matlab-*-

classdef FileUtils
    methods(Static)
        function yuv = convert420To444(yuv420)
            %takes in a vector of length R*C*1.5 and then returns the Y, U, V 444
            %format
            [rows cols] = size(yuv420{1});
            d1 = (yuv420{1}');
            d2 = (yuv420{2}');
            d3 = (yuv420{3}');
            data = [d1(:)' d2(:)' d3(:)'];
            yuv = zeros(rows,cols,3);
            yuv(:,:,1) = reshape(data(1:rows*cols),cols,rows)';
            newshape = reshape(data(rows*cols+1:rows*cols+1+(rows*cols/4)-1),cols/2,rows/2)';
            yuv(:,:,2) = double(imresize(newshape,2,'nearest'));
            newshape = reshape(data(rows*cols+(rows*cols/4)+1:end),cols/2,rows/2)';
            yuv(:,:,3) = double(imresize(newshape,2,'nearest'));
        end % convert420To444
        function yuv420 = convert444To420(yuv444)
            yuv420 = cell(1, 3);
            yuv420{1} = yuv444(:, :, 1);
            yuv420{2} = imresize(yuv444(:, :, 2), 0.5);
            yuv420{3} = imresize(yuv444(:, :, 3), 0.5);
        end % convert444To420
        function info = get_yuv420_info(filename)
            fid = fopen(filename);
            if(fid > 0)
                info.size  = FileUtils.GetSize(filename, 'x');
                if sum(info.size) == 0
                    % Try underscore if x doesn't work
                    info.size  = FileUtils.GetSize(filename, '_');
                end
                yuv_frame_size = 1.5*info.size(1)*info.size(2);
                fseek(fid, 0,'eof');
                info.frames = floor(ftell(fid)/yuv_frame_size);
                fclose(fid);
            else
                error('file %s can not be opened',filename);
            end
        end % get_yuv420_info
        function size = GetSize(filename, divide)
            size = [0 0];
            size_exp = regexp(filename,strcat('\d*',divide,'\d*'),'match');
            size_parse = regexp(size_exp, divide, 'split');
            if numel(size_parse) > 0
                rows = size_parse{end}(2);
                cols = size_parse{end}(1);
                size = [str2num(rows{1}),str2num(cols{1})];
            end
        end % GetSize
        function [width, height] = get_video_dims(filename)
            e = regexp(filename, '\.(\d+)(_|x)(\d+)\.', 'tokens');
            width = str2num(e{1}{1});
            height = str2num(e{1}{3});
        end

        function [yuv420] = readYUV420(filename,frame,depth)
            info = FileUtils.get_yuv420_info(filename);
            fid = fopen(filename);
            fseek(fid,1.5*prod(info.size)*(frame-1),'bof');
            if depth == 32
                data = fread(fid,1.5*prod(info.size),'float');
            else
                data = fread(fid,1.5*prod(info.size),'uchar');
            end
            fclose(fid);
            cols = info.size(2);
            rows = info.size(1);
            yuv420{1} = reshape(data(1:cols*rows), [cols, rows])';
            yuv420{2} = reshape(data((cols*rows+1):(cols*rows+((cols/2)*(rows/2)))), [cols/2, rows/2])';
            yuv420{3} = reshape(data(((cols*rows+((cols/2)*(rows/2)))+1):end), [cols/2,rows/2])';
        end % readYUV420
        function write420( fid, yuv420 )
            d1 = uint8(yuv420{1}');
            d2 = uint8(yuv420{2}');
            d3 = uint8(yuv420{3}');
            fseek(fid, 0, 'eof');
            fwrite(fid, d1(:), 'uint8');
            fwrite(fid, d2(:), 'uint8');
            fwrite(fid, d3(:), 'uint8');
        end % write420
        function convert420ToRGBRaw(filename_in, filename_out)
            file_info = get_yuv420_info(filename_in);
            fid = fopen(filename_out, 'a');
            if (fid < 0)
              error(sprintf('Failed to create file: %s', filename_out));
            end
            for frame = 1:file_info.frames
              yuv = FileUtils.readYUV420(filename_in, frame, 8);
              yuv = FileUtils.convert420To444(yuv);
              rgb = ycbcr2rgb(double(yuv)./255);
              d1 = zeros(1, prod(size(yuv)));
              d2(:,:,1) = uint8(255*rgb(:,:,1)');
              d3(:,:,1) = uint8(255*rgb(:,:,2)');
              d4(:,:,1) = uint8(255*rgb(:,:,3)');
              d1(1:3:end) = d2(:);
              d1(2:3:end) = d3(:);
              d1(3:3:end) = d4(:);
              fwrite(fid, d1(:), 'uint8');
            end
            fclose(fid);
        end % convert 8bit YUV 420 to 8bit RGB RAW
    end % methods
end % classdef
