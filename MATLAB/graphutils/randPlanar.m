function [edges, xy] = randPlanar(n)
%randPlanar(n) generates a random planar graph with n vertices
%
% this generates a planar graph using Delaunay triangulation
% the planar graph is not necessarily uniformly picked 
%
% Usage:   [edges, xy] = randPlanar(n)

xy = rand(n,2)*sqrt(n);
triangles = delaunay(xy(:,1), xy(:,2));

edges = [triangles(:,1:2); triangles(:,2:3); triangles(:,[1,3])];
edges = sortrows(edges);
edges = unique(edges,'rows');
