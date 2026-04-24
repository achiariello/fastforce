function [normals, centroids, areas] = computeFaceProperties(nodes, faces)
    % Compute normals, centroids, and areas of hexahedron faces
    % nodes: 3x8 matrix of node coordinates
    % faces: 6x4 matrix defining faces in terms of node indices
    % normals: 3x6 matrix of face normals
    % centroids: 3x6 matrix of face centroids
    % areas: 1x6 vector of face areas
    
    normals = zeros(3, size(faces, 1)); % Initialize normals
    centroids = zeros(3, size(faces, 1)); % Initialize centroids
    areas = zeros(1, size(faces, 1)); % Initialize areas

    for f = 1:size(faces, 1)
        % Nodes of current face
        face_nodes = nodes(:, faces(f, :));
        
        % Edge vectors on the face
        v1 = face_nodes(:, 2) - face_nodes(:, 1);
        v2 = face_nodes(:, 3) - face_nodes(:, 1);

        % Normal via cross product
        normal = cross(v1, v2);

        % Normalize normal vector
        normals(:, f) = normal / norm(normal);

        % Face centroid
        centroids(:, f) = mean(face_nodes, 2);

        % Face area
        areas(f) = norm(normal);
    end
end