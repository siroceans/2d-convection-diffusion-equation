%% MAE 623: CFD I
% Project 02
% Jorge Luis Mares Zamora
% Due 11/04/2025
clear 
clc
close all

%% Input parameters to run code. 

theta = 150; % [deg]
alpha = 1 * 10^(-6); % difussivity
Nx = 11; % Number of grid points in x
Ny = 11; % Number of grid points in ya
scheme = 1; % 1 for UD1, 2 for UD2, 3 for CDS

phi = convection_diffusionSolver(theta, alpha, Nx, Ny, scheme); 
phi = flipud(phi);

%% Beyond this point, the inputs are fixed to produce the results asked for in the report!
%% Discretization 1: UD1
phi_ud1_coarse = convection_diffusionSolver(theta, alpha, 11, 11, 1);
phi_ud1_medium = convection_diffusionSolver(theta, alpha, 21, 21, 1); 
phi_ud1_fine = convection_diffusionSolver(theta, alpha, 41, 41, 1); 

%% Discretization 2: UD2
phi_ud2_coarse = convection_diffusionSolver(theta, alpha, 11, 11, 2); 
phi_ud2_medium = convection_diffusionSolver(theta, alpha, 21, 21, 2); 
phi_ud2_fine = convection_diffusionSolver(theta, alpha, 41, 41, 2); 

%% Discretization 3: CDS
phi_cds_coarse = convection_diffusionSolver(theta, alpha, 11, 11, 3);
phi_cds_medium = convection_diffusionSolver(theta, alpha, 21, 21, 3); 
phi_cds_fine = convection_diffusionSolver(theta, alpha, 41, 41, 3);  

%% Plotting Results
plotReportResults(phi_ud1_coarse, phi_ud2_coarse, phi_cds_coarse, 11, 11)
plotReportResults(phi_ud1_medium, phi_ud2_medium, phi_cds_medium, 21, 21)
plotReportResults(phi_ud1_fine, phi_ud2_fine, phi_cds_fine, 41, 41)

plotResults3d(phi_ud1_fine, 41, 41, 1)
plotResults3d(phi_ud2_fine, 41, 41, 2)
plotResults3d(phi_cds_fine, 41, 41, 3)

%% Solver function
function [phi, iterations] = convection_diffusionSolver(theta, alpha, Nx, Ny, discretization)
    % check constraint on theta
    if (theta < 90) || (theta > 180)
        error('Theta must be in the following range: 90 < theta < 180')
    end

    % u and v are constant in the domain but dependent on theta
    u = cosd(theta);
    v = sind(theta);

    % initialize variables 
    dx = 1 / (Nx - 1); 
    dy = 1 / (Ny - 1); 
    tol = 10^(-10); % Error tolerance
    iteration = 0; % iteration counter

    % Direct method for the case of CDS discretization
    if discretization == 1 || discretization == 2
        gauss = true; 
        phi = zeros(Ny, Nx); % always index phi(j, i) <-> phi(row, column)
    else 
        gauss = false; % making sure not to go into the gauss seidel loop!
        N = Nx * Ny; 
        phi = zeros(N, 1); 
        
        % Coefficients for diagonal matrix
        % order - aw, as, ap, an, ae
        ap = (2*alpha/dx^2 + 2*alpha/dy^2); 
        ae = (u/(2*dx) - alpha/dx^2); 
        aw = (-u/(2*dx) - alpha/dx^2); 
        an = (v/(2*dy) - alpha/dy^2); 
        as = (-v/(2*dy) - alpha/dy^2); 

        diagonalArray = zeros(N, 5); 
        diagonalArray(:,1) = aw; 
        diagonalArray(:,2) = as; 
        diagonalArray(:,3) = ap; 
        diagonalArray(:,4) = an; 
        diagonalArray(:,5) = ae; 

        diagonalPositions = [-Ny, -1, 0, 1, Ny]; 
        A = spdiags(diagonalArray, diagonalPositions, N, N); 
        
        % Modifying the rows in the A matrix for the boundary conditions 
        left = 1:Ny; 
        right = Ny^2 - Ny + 1: Ny^2; 
        bottom = (0:Nx-1) * Nx + 1; 
        top = (1:Nx) * Nx; 

        for p = 1 : length(phi)
            if ismember(p, left)
                A(p, :) = 0; % resetting all coefficients to 0; 
                A(p, p) = 1; 
                A(p, p + Ny) = -2; 
                A(p, p + 2*Ny) = 1;  % p + 2 * Ny
            elseif ismember(p, right) % complete
                A(p, :) = 0; 
                A(p, p) = 1; 
            elseif ismember(p, top)
                A(p, :) = 0; 
                A(p, p) = 1; 
                A(p, p - 1) = -2; 
                A(p, p - 2) = 1; 
            elseif ismember(p, bottom) % complete
                A(p, :) = 0; 
                A(p, p) = 1;
            end
        end

        % Creating b vector
        b = zeros(size(phi)); 
        righttop = []; 
        for j = 0:Ny; 
            l = (Nx - 1) * Ny + j; 
            if (j * dy) > 0.25
               righttop = [righttop, l]; 
            end
        end

        for p = 1:length(b)
            if ismember(p, righttop)
                b(p) = 1; 
            end
        end
        phi = A \ b; 
        phi = reshape(phi, [Nx, Ny]); 
    end
    
    % Using gauss seidel
    while gauss
        phi_old = phi; % storing previous iterate for error calculation!
        iteration = iteration + 1;

        % Updating values at interior nodes
        for i = 2:(Nx - 1)
            for j = 2:(Ny - 1)
                % UD1
                if discretization == 1 % using UD1
                    numerator = -(u/dx - alpha/dx^2)*phi(j,i+1) - (-v/dy - alpha/dy^2)*phi(j-1,i) - (-alpha/dx^2)*phi(j,i-1) - (-alpha/dy^2)*phi(j+1,i);
                    denominator = (v/dy - u/dx + 2*alpha/dx^2 + 2*alpha/dy^2); 
                    phi(j, i) = numerator / denominator; 
                % UD2
                elseif discretization == 2 % using UD2
                    % if there is no room to index j - 2 or i + 2 
                    if j == 2 || i == (Nx-1)
                        % we do a UD1 discretization 
                        numerator = -(u/dx - alpha/dx^2)*phi(j,i+1) - (-v/dy - alpha/dy^2)*phi(j-1,i) - (-alpha/dx^2)*phi(j,i-1) - (-alpha/dy^2)*phi(j+1,i);
                        denominator = (v/dy - u/dx + 2*alpha/dx^2 + 2*alpha/dy^2); 
                        phi(j, i) = numerator / denominator; 
                    else % if not on those edge cases, we can safely do UD2!!!
                        numerator = (2*u/dx - alpha/dx^2)*phi(j, i+1) + (-u/(2*dx))*phi(j,i+2) + (v/(2*dy))*phi(j-2,i) + (-2*v/dy - alpha/dy^2)*phi(j-1,i) + (-alpha/dx^2)*phi(j,i-1) + (-alpha/dy^2)*phi(j+1,i); 
                        denominator = 3*u/(2*dx) - 3*v/(2*dy) - 2*alpha/dx^2 - 2*alpha/dy^2; 
                        phi(j, i) = numerator / denominator;
                    end

                % CDS -- cannot be solved using gauss seidel, solved using a direct method instead!
                % see updated code above
                %elseif discretization == 3
%                    numerator = (u/(2*dx) - alpha/dx^2)*phi(j,i+1) + (-u/(2*dx) - alpha/dx^2)*phi(j,i-1) + (v/(2*dy) - alpha/dy^2)*phi(j+1,i) + (-v/(2*dy) - alpha/dy^2)*phi(j-1, i); 
                    %denominator = -(2*alpha/dx^2 + 2*alpha/dy^2); 
                    %phi(j, i) = numerator / denominator;
                end
            end
        end

        % Updating boundaries!
        for i = 1 : Nx 
            % SOUTH boundary
            phi(1, i) = 0; 
            % NORTH boundary
            phi(Ny, i) = 2*phi(Ny-1, i) - phi(Ny-2, i); % updated for linear extrapolation
        end

        for j = 1:Ny 
            % EAST boundary
            if (j * dy) < 0.25 
                phi(j, Nx) = 0; 
            else 
                phi(j, Nx) = 1; 
            end

            % WEST boundary
            phi(j, 1) = 2 * phi(j, 2) - phi(j, 3); % updated for linear extrapolation
        end

        % Convergence criteria
        e = abs((phi - phi_old)); % Magnitude approximate error! 
        if max(max(e)) < tol %| iteration > 100 % limiting the number of iterations!
            break
        end
    end
end


%% Helper functions
function plotResults3d(phi, Nx, Ny, no)
    %phi = flipud(phi); % change from matlab upside down ordering!
    x = linspace(0, 1, Nx); 
    y = linspace(0, 1, Ny); 
    [X,Y] = meshgrid(x, y); 

    % title name conditional
    if no == 1
        scheme = "UD1"; 
    elseif no == 2
        scheme = "UD2"; 
    elseif no == 3
        scheme = "CDS"; 
    end

    figure()
    surf(X, Y, phi)
    title(['2D Convection/Difussion Equation Using' scheme])
    xlabel('x direction')
    ylabel('y direction')
    zlabel('z direction')
end

function plotReportResults(phi1, phi2, phi3, nx, ny)
    phi1 = flipud(phi1); 
    phi2 = flipud(phi2); 
    phi3 = flipud(phi3); 
    
    figure()
    yvalues = 0:(1/(ny-1)):1;
    phivals = phi1(:, (nx-1)/2)';
    phivals2 = phi2(:, (nx-1)/2)'; 
    phivals3 = phi3(:, (nx-1)/2)'; 
    width = 2; 
    plot(yvalues, phivals,'k', 'LineWidth', width)
    hold on
    plot(yvalues, phivals2, 'b', 'LineWidth', width)
    plot(yvalues, phivals3, 'm', 'LineWidth', width)
    hold off
    
    title(['Phi vs. y at x = 0.5 for an ' num2str(nx) ' x ' num2str(ny) 'grid resolution'])
    xlabel 'y values' 
    ylabel 'Phi values'
    grid on 
    legend('UD1 discretization of the convection term.', 'UD2 discretization of the convection term.', 'CDS discretization of the convection term.', 'location', 'best')
end 
