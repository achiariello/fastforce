clear
close all
% addpath('..\CodiciFL\ToAndrea\hexahedron')

% Constants
mu0 = 4 * pi * 1e-7; % Permeability of free space
Tol = 1e-6; % Convergence tolerance
maxit = 100; % Maximum iterations
NumberDipoles = 3; % Number of dipoles
delta = 1e-6; % Finite difference step

% Node indices for each face
faces = [
    1, 4, 3, 2; % Face 1
    5, 6, 7, 8; % Face 2
    1, 2, 6, 5; % Face 3
    2, 3, 7, 6; % Face 4
    3, 4, 8, 7; % Face 5
    4, 1, 5, 8  % Face 6
];

T_faces = zeros(3, 3, 6, NumberDipoles); % Maxwell tensor for each face and dipole
B_faces_store = zeros(3, 6, NumberDipoles); % Field for each face and dipole
B_totals_store = zeros(3, NumberDipoles); % Field at each dipole barycenter

% Input data
P1 = [12.861, -0.191+0.152, -0.3545+0.271];
P2 = [12.868, -0.0815+0.152, -0.145+0.271];
P3 = [13.115, -0.181+0.152, -0.317+0.271];

% Circular coil source along the axis
curr_spira=1.7e6;
x0=11;
R=7;

B_source = @(x) (mu0 * curr_spira * R^2) ./ (2 * ((x-x0).^2 + R^2).^(3/2));

dB_source_dx=@(x) - (3 * mu0* curr_spira * R^2 .* (x - x0)) ./ (2 * ((x - x0).^2 + R^2).^(5/2));

dBsource_dx_store= zeros(1,NumberDipoles);

% Uniform magnetic field in the fixator
B1 = [B_source(P1(1)),0, 0]; 
B2 = [B_source(P2(1)),0, 0]; 
B3 = [B_source(P3(1)),0, 0]; 

B0=[B1; B2; B3];

fNL= @(B)  (1/(4*pi*1e-7))*(1-1/MB_fixator_ur(B))*B;

% Initialization
M = 0.*rand(3, NumberDipoles); % Initial arbitrary magnetizations
Error = 1e6; % Initial error
it = 0; % Iteration counter
Error_vet = []; % Error at each iteration

% Dipole positions
positions = [P1; P2; P3];

% Hexahedron dimensions
dLato_all = [120 178 75;
    150 95 220;
    95 220 150]*1e-3;

% Iterative solution
while (it <= maxit && Error > Tol)
    it = it + 1;
    M_new = zeros(3, NumberDipoles); % Updated magnetizations
    nablaB_other = zeros(3, 3, NumberDipoles); % Gradient matrix

    for i = 1:NumberDipoles
        % Field due to all other dipoles
        B_other = [0; 0; 0];
        
        for j = 1:NumberDipoles
            dLatoj=dLato_all(j,:);
            if j ~= i
                baricenter_j = positions(j, :)'; % Barycenter of hexahedron j
                nodes_j = generateHexNodes(baricenter_j, dLatoj); % Nodes of hexahedron j

                P_i = positions(i, :)'; % Field point
                J = [0; 0; 0]; % Current
                M_j = M(:, j); % Magnetization of hexahedron j

                [~, ~, ~, B_j] = hexahedron(nodes_j, P_i, J, M_j);

                B_other = B_other + B_j;

                % Gradient of B_other by finite differences
                for dim = 1:3
                    P_shift_p = P_i;
                    P_shift_m = P_i;
                    P_shift_p(dim) = P_shift_p(dim) + delta;
                    P_shift_m(dim) = P_shift_m(dim) - delta;

                    [~, ~, ~, B_j_shifted_p] = hexahedron(nodes_j, P_shift_p, J, M_j);
                    [~, ~, ~, B_j_shifted_m] = hexahedron(nodes_j, P_shift_m, J, M_j);

                    dB_j = (B_j_shifted_p - B_j_shifted_m) /(2* delta);
                    nablaB_other(:, dim, i) = nablaB_other(:, dim, i) + dB_j;
                end
            end
        end

        baricenter_i = positions(i, :)'; % Barycenter of hexahedron i
        dLatoi=dLato_all(i,:);
        nodes_i = generateHexNodes(baricenter_i, dLatoi); % Nodes of hexahedron i

        P_i = positions(i, :)'; % Field point
        dBsource_dx_store(:,i)=dB_source_dx(P_i(1));

        M_i = M(:, i); % Magnetization of hexahedron i
        [~, ~, ~, B_self] = hexahedron(nodes_i, P_i, J, M_i);

        B_total = B0(i,:)' + B_other + B_self;
        B_totals_store(:,i)=B_total;

        M_new(:, i) = fNL(B_total);
    end

    % Error computation
    Error = 0;
    for i = 1:NumberDipoles
        Error = Error + norm(M(:, i) - M_new(:, i)) / (norm(M_new(:, i)) + eps);
    end
    Error_vet(it) = Error;

    M = M_new;
end

% Display results
fprintf('Iterations: %d\n', it);
fprintf('Final Error: %.5e\n', Error);
disp('Final Magnetizations (M1 and M2):');
disp(M);

% Plot error convergence
figure;
semilogy(1:it, Error_vet, '-o');
xlabel('Iteration');
ylabel('Error (log scale)');
title('Convergence of Error');
grid on;

% Resulting forces on the hexahedra
forces = zeros(3, NumberDipoles); % Total forces on each hexahedron
areas_store=zeros(6, NumberDipoles);

for i = 1:NumberDipoles
    baricenter_i = positions(i, :)'; % Barycenter of dipole i
    dLatoi=dLato_all(i,:);
    nodes_i = generateHexNodes(baricenter_i, dLatoi); % Nodes of hexahedron i
    
    forces(:, i) = computeForces(nodes_i, faces, T_faces(:, :, :, i));

    % Face normals and centroids
    [normals, centroids,areas] = computeFaceProperties(nodes_i, faces);
    areas_store(:,i)=areas';

    figure;
    hold on;
    scatter3(nodes_i(1, :), nodes_i(2, :), nodes_i(3, :), 'filled', 'r');

    for f = 1:6
        face_nodes = nodes_i(:, faces(f, :));
        fill3(face_nodes(1, :), face_nodes(2, :), face_nodes(3, :), 'g', 'FaceAlpha', 0.2);
    end

    quiver3(centroids(1, :), centroids(2, :), centroids(3, :), ...
            normals(1, :), normals(2, :), normals(3, :), 0.1, 'k', 'LineWidth', 2);

    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    axis equal;
    grid on;
    title(['Hexahedron ', num2str(i), ' with Face Normals']);
    hold off;
end

dV=dLato_all(1,1)*dLato_all(1,2)*dLato_all(1,3);
Fdipole=zeros(3, NumberDipoles);

for i = 1:NumberDipoles
    baricenter_i = positions(i, :)';
    for j = 1:NumberDipoles
        baricenter_j = positions(j, :)';
        if (i~=j)
        Fdipole(:,i) =Fdipole(:,i)+ dipole_force(M(:,i).*dV, M(:,j).*dV, (baricenter_i-baricenter_j));
        end
    end
end

disp('Dipole formula forces for each dipole:');
disp(Fdipole)

% Kelvin force initialization
Fkelvin = zeros(3, NumberDipoles);
Fkelvin_sources= zeros(3, NumberDipoles);

% Dipole volumes
dV_all = prod(dLato_all, 2);

for i = 1:NumberDipoles
    % H_ext gradient at the barycenter
    H_ext_grad = nablaB_other(:, :, i) /mu0;

    % Force contribution for dipole i
    Fkelvin(:, i) = mu0 * (M(:, i)' * H_ext_grad)' * dV_all(i);
    Fkelvin_sources(1,i)= mu0 * (M(1, i)' * dBsource_dx_store(i)/mu0)' * dV_all(i);
end

disp('Kelvin forces for each dipole without external source contribution:');
disp(Fkelvin);

disp('Kelvin forces for each dipole: external source contribution only:');
disp(Fkelvin_sources);

disp('Total forces')
disp((Fkelvin+Fkelvin_sources))