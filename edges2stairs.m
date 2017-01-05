function [X,Y] = edges2stairs(edges,counts)

X = reshape([reshape(edges,1,[]);reshape(edges,1,[])],[],1);
Y = reshape([reshape(counts,1,[]);reshape(counts,1,[])],[],1);
Y = [0;Y;0];