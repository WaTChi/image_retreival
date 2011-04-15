% transfer 4th floor cory images based on filtered list

fn = 'imagesUsed_7240131.texture';

src = '\\bam\Shared Folder\Data Collection\May 4 2010\20100504-3\';
dest = 'Z:\cory\';

files = textread(fn,'%s');
files = files(3:2:end);

for k=1:length(files)
    file = [files{k}(13:end),'.jpg'];
    copyfile([src,file],[dest,file]);
end