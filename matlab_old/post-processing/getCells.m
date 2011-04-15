function [cells,cell_lats,cell_lons] = getCells

% [cell_files,cell_lats,cell_lons] = getCells(query,vote_dir)
% 
% DESCRIPTION
%   This function reads through the cell map  and returns the names of
%   the cells as well as their latitude and longitude locations.

[~,cells] = textread('C:\matlab_local\cellmap.txt','%d%s');
comma_idx = strfind(cells,',');
ncells = length(cells);
cell_lats = zeros(ncells,1);
cell_lons = zeros(ncells,1);
for k=1:ncells
    cell_lats(k) = str2double( ...
        cells{k}( 1:comma_idx{k}-1 ) );
    cell_lons(k) = str2double( ...
        cells{k}( comma_idx{k}+1:end ) );
end