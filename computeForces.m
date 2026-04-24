function force_total = computeForces(nodes, faces, T_faces)
    % Compute total force acting on a hexahedron
    % nodes: 3x8 matrix of hexahedron nodes
    % faces: 6x4 matrix defining the faces
    % T_faces: Maxwell tensor for each face
    % force_total: 3x1 total force vector

    % Compute face properties
    [normals, ~, areas] = computeFaceProperties(nodes, faces);

    % Sum of forces over all faces
    force_total = [0; 0; 0];
    for f = 1:size(faces, 1)
        % Maxwell tensor on face f
        T_face = T_faces(:, :, f);

        % Face normal
        normal_f = normals(:, f);

        % Face area
        area_f = areas(f);

        % Force on face
        force_f = T_face * normal_f * area_f;

        % Accumulate total force
        force_total = force_total + force_f;
    end
end