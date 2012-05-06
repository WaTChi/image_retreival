function [cellgroups] = groupFuzzy(lats,lons,cell_lats,cell_lons,cell_dist)

% [cellgroups] = groupFuzzy(lats,lons,cell_lats,cell_lons,cell_dist)
% 
% DESCRIPTION
%   This function creates groupings of fuzzy points depending on which
%   cells the points search over based on the maximum cell search distance.
% 
% INPUTS
%   lats:       Fuzzy point latitude coordinates
%   lons:       Fuzzy point longitude coordinates
%   cell_lats:  Cell center latitude coordinates
%   cell_lons:  Cell center longitude coordinates
%   cell_dist:  Threshold distance between cell center and fuzzy point to
%               trigger search. If equal to 0, this is no combination and
%               we choose closest cell.
% 
% OUTPUTS
%   cellgroups: Structure array containing the cell grouping information
%       .idx:   Array indicating the different indices of the cells 
%               searched over in this grouping
%       .lats:  Array of fuzzy point latitude coordinates in this grouping
%       .lons:  Array of fuzzy point longitude coordinates in this grouping
%       .npts:  Number of fuzzy points in this grouping

F = length(lats);
C = length(cell_lats);

mat_lats = repmat( lats, [1,C] );
mat_lons = repmat( lons, [1,C] );
cell_lats = repmat( cell_lats', [F,1] );
cell_lons = repmat( cell_lons', [F,1] );

distances = latlonDistance( mat_lats,mat_lons , cell_lats,cell_lons );

% Get search flags based on combination method
if cell_dist == 0 % no combination
    min_distances = repmat( min( distances,[],2 ) , [1,C] );
    search_flags = min_distances == distances;
else % combination
    search_flags = distances < cell_dist;
end

k = 1;
while ~isempty(lats)
    cellgroups(k).idx = search_flags(1,:);
    group_idx = repmat( cellgroups(k).idx , [length(lats),1] );
    fuzzy_mask = logical( prod( double( ~xor(search_flags,group_idx) ) , 2 ) );
    cellgroups(k).lats = lats(fuzzy_mask);
    cellgroups(k).lons = lons(fuzzy_mask);
    cellgroups(k).npts = sum(fuzzy_mask);
    search_flags(fuzzy_mask,:) = [];
    lats(fuzzy_mask) = [];
    lons(fuzzy_mask) = [];
    if any(cellgroups(k).idx)
        k=k+1;
    end
end
cellgroups = cellgroups(1:k-1);