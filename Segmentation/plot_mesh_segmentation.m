function [seg,V,F] = plot_mesh_segmentation(meshFile, segFile)
% Given a mesh file (.off or .obj) and a segmentation
% (.seg) file, it plots the mesh with different colors
% for the different segments and outputs the number of 
% segments
%
% Example usage:
%  f = plot_mesh_segmentation('../MeshsegBenchmark-1.0/data/off/285.off','../MeshsegBenchmark-1.0/data/seg/Benchmark/285/285_0.seg')
%
% Input:
% meshFile: path to mesh file (.off or .obj)
% segFile: path to segmentation file (.seg)
%
% Output:
% seg: vector containing segment index for each face


fileID = fopen(segFile,'r');
formatSpec = '%d';
seg = fscanf(fileID,formatSpec);
fclose(fileID);
seg = seg+1;
max_label = max(seg);
colors = rand(max_label,3);

[V,F] = load_mesh(meshFile);
for k=1:max(seg)
    tsurf(F(seg==k,:),V,'FaceColor',colors(k,:))
    hold on
end
axis equal