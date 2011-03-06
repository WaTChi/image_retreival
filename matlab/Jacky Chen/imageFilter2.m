function [] = imageFilter2(dataset, varargin)

%   imageFilter(dataset, option) filters an entire image/feature set based
%   on the texture file, so that at the end, only those images/features
%   listed in the texture file remain in the folder.
%
%   Input:
%       - dataset: the folder in which the images/features reside. Needs to
%       be a string. e.g. '20100504_1' The directory for the code will then
%       change to .\20100504_1\.
%       - varargin: an additional argument of 'sift' will allow the code to
%       filter features. Otherwise, it defaults to filtering images.
%
%   Output:
%       - Nothing is returned. Only the images/features specified by the 
%       texture file would remain inside the folder.
%

format = '*.jpg';
if size(varargin, 2) == 1 && strcmpi(varargin{1}, 'sift')
    format = '*.feat';
    option = varargin{1};
end

num = 6;

dir = ['.\', dataset, '\']; % must end in file separator (e.g. '\'); this 
            % is an input directory that must consist of *.mat files 
            % specifying SIFT features according to a format that Matt 
            % Carlberg was using

siftFeatures = ls(fullfile(dir, format)); 
siftFeatures = [repmat(dir, size(siftFeatures, 1), 1), siftFeatures];
N = size(siftFeatures, 1); % N = total number of documents/images

fileName = [dir, ls(fullfile(dir, '*.texture'))];
fid = fopen(fileName, 'r');
textscan(fid, '%d', 1);
images = textscan(fid, '%*d%*s%s', 'delimiter', '/');
images = images{1};
K = size(images, 1);

k = 1;
chosen = images{k};
e = num2str(chosen(end-num+1:end));
for i = 1:N
    path = siftFeatures(i,:);
    if strcmpi(option, 'sift')
        d = textscan(path(end-num+1-5:end), '%d', 'delimiter', '.');
    else
        d = textscan(path(end-num+1-4:end-4), '%d', 'delimiter', '.');
    end
    d = d{1};
    if isempty(strfind(path, chosen)) % image not needed; can be thrown
        while (d > e) % advance to next wanted
            k = k + 1;
            if k > K % all images in texture file seen; delete the rest
                for j = i:N
                    delete(siftFeatures(j,:))
                end
                return
            end
            chosen = images{k};
            e = num2str(chosen(end-num+1:end));
        end
        if e > d
            delete(path)
        end
    else % not thrown away and move on
        k = k + 1;
        if k > K % all images in texture file seen; delete the rest
            for j = i+1:N
                delete(siftFeatures(j,:)) 
            end
            return
        end
        chosen = images{k};
        e = num2str(chosen(end-num+1:end));
    end
end
fclose(fid);