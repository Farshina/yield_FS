function [idx,scores] = findmrmr(X,y)
% perform MRMR in the dataset
[idx,scores] = fsrmrmr(X,y);

% normalize scores
scores = normalize(scores, 'range');
end
