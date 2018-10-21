function plotCoordinatesTransform(T, scale)
% PlotCoordinate: Plot coordinates in the XYZ-RGB sequence
origin = T(1:3,4);
unit_vecs = T(1:3, 1:3);
axis_x = unit_vecs(:,1);
axis_y = unit_vecs(:,2);
axis_z = unit_vecs(:,3);
x = origin(1); y = origin(2); z = origin(3);
quiver3(x,y,z,axis_x(1),axis_x(2),axis_x(3),scale,'color','r'); hold on
quiver3(x,y,z,axis_y(1),axis_y(2),axis_y(3),scale,'color','g')
quiver3(x,y,z,axis_z(1),axis_z(2),axis_z(3),scale,'color','b')
end

