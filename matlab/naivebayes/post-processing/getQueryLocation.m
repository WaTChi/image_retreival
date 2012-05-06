function [qLat,qLon] = getQueryLocation(query,set)

if strcmp(set,'query1') || strcmp(set,'query2') || strcmp(set,'query3')
    
    % parse file name format
    idx = strfind(query,',');
    qLat = str2double(query(idx(1)+1:idx(2)-1));
    qLon = str2double(query(idx(2)+1:end-8));
    
else % anything but sets 1,2,3
    
    % read xml contents
    qdir = ['Z:\',set,'\'];
    fn_base = query(1:end-8);
    fn = [qdir,'imageotag_',fn_base,'.xml'];
    fid = fopen(fn);
    if fid<0
        fclose(fid);
        error(['File not found: ',fn])
    end
    text = fscanf(fid,'%s');
    fclose(fid);
    
    % search for latitude
    idx1 = strfind(text,'<gps_latitude>');
    idx2 = strfind(text,'</gps_latitude>');
    qLat = str2double(text(idx1+14:idx2-1));
    
    % search for longitude
    idx1 = strfind(text,'<gps_longitude>');
    idx2 = strfind(text,'</gps_longitude>');
    qLon = str2double(text(idx1+15:idx2-1));
    
end
