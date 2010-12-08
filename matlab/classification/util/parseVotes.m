function [qVote,qPhoto] = parseVotes(filename,vote_dir)

if nargin<2
    vote_dir = 'E:\q3results\';
end

vote_data = textread([vote_dir,filename],'%s');
vote_data = reshape(vote_data,[2,length(vote_data)/2])';
qVote = str2num(strvcat(vote_data(:,1)));
qPhoto = vote_data(:,2);
for k=1:length(qPhoto)
    qPhoto{k} = qPhoto{k}(1:end-8);
end

end