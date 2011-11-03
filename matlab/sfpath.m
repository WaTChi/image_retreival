inputdir = 'Z:\ah\sfpath\';
outfile = 'Z:\ah\sfpath.txt';

inputfiles = struct2cell(dir(inputdir));
inputfiles = inputfiles(1,:)';
geotag_idx = ~cellfun('isempty',strfind(inputfiles,'imageotag_'));
inputfiles = inputfiles(geotag_idx);

fout = fopen(outfile,'w');

for k=1:length(inputfiles)
    
    fin = fopen([inputdir,inputfiles{k}]);
    text = fscanf(fin,'%s');
    fclose(fin);

    % search for latitude
    idx1 = strfind(text,'<gps_latitude>');
    idx2 = strfind(text,'</gps_latitude>');
    lat = str2double(text(idx1+14:idx2-1));

    % search for longitude
    idx1 = strfind(text,'<gps_longitude>');
    idx2 = strfind(text,'</gps_longitude>');
    lon = str2double(text(idx1+15:idx2-1));
    
    fprintf(fout,'%d, ,%.9f,%.9f\n',k,lon,lat);
    
end

fclose(fout);