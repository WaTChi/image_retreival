path ='~/Desktop/query/project/src/tutorial/mall_query'
files = dir('~/Desktop/query/project/src/tutorial/mall_query/*.jpg');
for idx = 1:numel(files)
    file = files(idx);
    disp(file.name);
    correct_image(path, file.name);
end