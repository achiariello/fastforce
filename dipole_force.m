function F = dipole_force(m1, m2, r)
    % Compute force between two magnetic dipoles
    % m1, m2: magnetic dipole moments [Mx, My, Mz]
    % r: relative position vector [rx, ry, rz]
    % F: force vector [Fx, Fy, Fz]

    % Physical constant
    mu0 = 4 * pi * 1e-7; % Vacuum permeability (H/m)

    % Distance and unit vector
    r_mag = norm(r); % Magnitude of r
    r_hat = r / r_mag; % Unit vector of r

    % Dot products
    m1_dot_rhat = dot(m1, r_hat);
    m2_dot_rhat = dot(m2, r_hat);
    m1_dot_m2 = dot(m1, m2);

    % Force terms
    term1 = m2_dot_rhat * m1;
    term2 = m1_dot_rhat * m2;
    term3 = m1_dot_m2 * r_hat;
    term4 = -5 * m1_dot_rhat * m2_dot_rhat * r_hat;

    % Force computation
    F = (3 * mu0 / (4 * pi * r_mag^4)) * (term1 + term2 + term3 + term4);
end