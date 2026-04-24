function nodes = generateHexNodes(baricenter, dLato)
    % Generate nodes of a hexahedron centered at the given barycenter
    % baricenter: [x, y, z] position of the barycenter
    % dLato: side lengths of the hexahedron [Lx, Ly, Lz]
    % nodes: 3x8 matrix containing node coordinates
    
    % Half side lengths
    halfLx = dLato(1) / 2; 
    halfLy = dLato(2) / 2; 
    halfLz = dLato(3) / 2;

    % Cube vertices centered at the origin
    baseCube = [
        -halfLx,  halfLx,  halfLx, -halfLx, -halfLx,  halfLx,  halfLx, -halfLx;
        -halfLy, -halfLy,  halfLy,  halfLy, -halfLy, -halfLy,  halfLy,  halfLy;
        -halfLz, -halfLz, -halfLz, -halfLz,  halfLz,  halfLz,  halfLz,  halfLz
    ];

    % Translation to the desired barycenter
    nodes = baseCube + baricenter(:);
end