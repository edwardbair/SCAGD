function [x11,y11,dx,dy,inProj] = cornerCoords(R)
% corner coordinates from raster reference
if contains(class(R),'MapCellsReference')
    inProj = true;
    [x11,y11] = intrinsicToWorld(R,1,1);
    if strcmp(R.ColumnsStartFrom,'north')
        dy = -R.CellExtentInWorldY;
    else
        dy = R.CellExtentInWorldY;
    end
    if strcmp(R.RowsStartFrom,'west')
        dx = R.CellExtentInWorldX;
    else
        dx = -R.CellExtentInWorldX;
    end
elseif contains(class(R),'GeographicCellsReference')
    inProj = false;
    [y11,x11] = intrinsicToGeographic(R,1,1);
    if strcmp(R.ColumnsStartFrom,'north')
        dy = -R.CellExtentInLatitude;
    else
        dy = R.CellExtentInLatitude;
    end
    if strcmp(R.RowsStartFrom,'west')
        dx = R.CellExtentInLongitude;
    else
        dx = -R.CellExtentInLongitude;
    end
else
    error('raster reference class %s unrecognized',class(R));
end
end