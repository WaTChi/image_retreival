function [photos,votes] = combineCells(vote_dir,comb_dir)

earth_radius = 6371000; % meters
comb_files = dir(comb_dir);
comb_files = strvcat(comb_files.name);
comb_nums = unique(str2double(cellstr(comb_files(3:end,5:8))));
comb_files = cellstr(comb_files(3:end,:));
query_files = dir(vote_dir);
query_files = strvcat(query_files.name);
query_files = query_files(3:end,:);
query_nums = unique(str2double(cellstr(query_files(:,5:8))));
% tmp = [comb_nums;8937]-query_nums;
% tmp2 = find(tmp);
% tmp2 = tmp2(1)
% query_nums(29)
% query_nums(30)
% comb_nums(29)
% comb_nums(30)
% query_nums(1)
% comb_files{1}
% query_nums(84)
% comb_files{84}
photos = cell(0,1);
votes = zeros(0,1);
for k=1:length(query_nums)
    disp(num2str(query_nums(k)))
    file_idx = strmatch(num2str(query_nums(k)),query_files(:,5:8));
    files = cellstr(query_files(file_idx,:));
    coords = files{1};
    comma_idx = strfind(coords,',');
    qLat = pi/180*str2double(coords(comma_idx(1)+1:comma_idx(2)-1));
    qLong = pi/180*str2double(coords(comma_idx(2)+1:comma_idx(3)-9));
    ncells = length(files);
%     range_idx = [];
    for j=1:ncells
        start_idx = strfind(files{j},',');
        start_idx = 1+start_idx(end);
        end_idx = -1+strfind(files{j},'.res');
        cell_dist = str2double(files{j}(start_idx:end_idx));
        if cell_dist < 336.6
%             range_idx = [range_idx,j];
            [qVote,qPhoto] = parseVotes(files{j},vote_dir);
            photos = [photos;qPhoto];
            votes = [votes;qVote];
        end
    end
    images = cell(0,1);
    number = zeros(0,1);
    while ~isempty(photos)
        images = [images;photos(1)];
        idx = strcmp(photos,photos(1));
        number = [number;sum(votes(idx))];
        photos(idx) = [];
        votes(idx) = [];
    end
    [number,sort_idx] = sort(number,'descend');
    images = images(sort_idx);
    fid = fopen([comb_dir,comb_files{k}],'w');
    for j=1:length(images)
        comma_idx = strfind(images{j},',');
        dbLat = pi/180*str2double(images{j}(1:comma_idx-1));
        dbLong = pi/180*str2double(images{j}(1+comma_idx:end-13));
        dbDist = earth_radius * acos(sin(qLat)*sin(dbLat) + ...
                    cos(qLat)*cos(dbLat)*cos(qLong-dbLong) );
        if dbDist < 100;
            fprintf(fid,'%d\t%s\n',number(j),images{j});
        end
    end
    fclose(fid);
end
%     files = files(range_idx);
%     start_idx = start_idx(:,end)+1;
%     end_idx = cell2mat(strfind(files,'.res'))-1;
    