%% Optimised Magnetic Field Analysis for Different Coil Configurations - Version 2.0
% Author: Adil Wahab Bhatti
% Version: 2.0
% Date: 2025-06-10
% Description: Advanced vectorised analysis of magnetic fields from Single, 
%             Helmholtz, Maxwell, and Double-Pair Helmholtz coil configurations
% File: magnetic_field_analyser_v2.m

clear; clc; close all;

%% Constants and Parameters
mu0 = 4 * pi * 1e-7;   % Permeability of free space (T*m/A)
N = 1300;              % Number of turns in each coil
I = 0.55;              % Current in Amperes
r = 0.1409;            % Radius of the coils in meters
z = linspace(-r, r, 1000); % z-axis range

%% Coil Configuration Structure
% Define all coil configurations in a structured way
configs = struct();

% Single Coil
configs.single = struct(...
    'name', 'Single Coil', ...
    'positions', [0], ...
    'radii', [r], ...
    'turns', [N], ...
    'color', '#0072BD');

% Helmholtz Coils
theta_H = acos(1/sqrt(5));
configs.helmholtz = struct(...
    'name', 'Single-Pair Helmholtz Coil', ...
    'positions', [-r*cos(theta_H), r*cos(theta_H)], ...
    'radii', [r*sin(theta_H), r*sin(theta_H)], ...
    'turns', [N, N], ...
    'color', '#D95319');

% Maxwell Coils
configs.maxwell = struct(...
    'name', 'Maxwell Coil', ...
    'positions', [-r*sqrt(3/7), 0, r*sqrt(3/7)], ...
    'radii', [r*sqrt(4/7), r, r*sqrt(4/7)], ...
    'turns', [floor(N*49/64), floor(N), floor(N*49/64)], ...
    'color', '#77AC30');

% Double-Pair Helmholtz Coils
theta_H1 = 40.09;
theta_H2 = 73.43;
configs.double_helmholtz = struct(...
    'name', 'Double-Pair Helmholtz Coil', ...
    'positions', [-r*cosd(theta_H1), -r*cosd(theta_H2), r*cosd(theta_H2), r*cosd(theta_H1)], ...
    'radii', [r*sind(theta_H1), r*sind(theta_H2), r*sind(theta_H2), r*sind(theta_H1)], ...
    'turns', [floor(N*0.68212), floor(N), floor(N), floor(N*0.68212)], ...
    'color', '#7E2F8E');

%% Core Function: Calculate Magnetic Field
function [B_total, B_individual] = calculate_magnetic_field(z, positions, radii, turns, mu0, I)
    % Vectorised calculation of magnetic field from multiple coils
    % Returns: B_total - combined magnetic field, B_individual - field from each coil
    
    B_total = zeros(size(z));
    B_individual = zeros(length(positions), length(z));
    
    for i = 1:length(positions)
        z_diff = z - positions(i);
        denominator = (radii(i)^2 + z_diff.^2).^(3/2);
        B_coil = (mu0 * turns(i) * I * radii(i)^2) ./ (2 * denominator);
        B_individual(i, :) = B_coil;
        B_total = B_total + B_coil;
    end
end

%% Calculate Fields for All Configurations
field_names = fieldnames(configs);
results = struct();

for i = 1:length(field_names)
    config_name = field_names{i};
    config = configs.(config_name);
    
    % Calculate magnetic field
    [B_total, B_individual] = calculate_magnetic_field(...
        z, config.positions, config.radii, config.turns, mu0, I);
    
    % Store results
    results.(config_name) = struct(...
        'B_total', B_total, ...
        'B_individual', B_individual, ...
        'B_max', max(B_total), ...
        'config', config);
end

%% Enhanced Plotting Function with Dual Format Output
function create_plot(z, r, results, config_name, plot_filename, show_decomposition, use_assigned_colors)
    if nargin < 6
        show_decomposition = false;
    end
    if nargin < 7
        use_assigned_colors = false;
    end
    
    % Create figure with invisible property initially
    fig = figure('Units', 'centimeters', 'Position', [0, 0, 8, 6], 'Visible', 'off');
    hold on;
    
    result = results.(config_name);
    config = result.config;
    
    % Get MATLAB default colors
    default_colors = get(gca, 'ColorOrder');
    
    if show_decomposition && size(result.B_individual, 1) > 1
        % Plot individual coil contributions
        for i = 1:size(result.B_individual, 1)
            if use_assigned_colors
                line_color = config.color;
            else
                % Use MATLAB default colors, cycling through them
                color_idx = mod(i-1, size(default_colors, 1)) + 1;
                line_color = default_colors(color_idx, :);
            end
            
            plot(z/r, result.B_individual(i,:)/result.B_max, '--', ...
                'Color', line_color, 'LineWidth', 0.5, ...
                'DisplayName', sprintf('Coil %d', i));
        end
    end
    
    % Plot combined field
    if use_assigned_colors
        line_color = config.color;
    elseif strcmp(config_name, 'single')
        line_color = 'k';  % Black for single coil
    else
        % Use the next color in sequence for the combined field
        color_idx = mod(size(result.B_individual, 1), size(default_colors, 1)) + 1;
        line_color = default_colors(color_idx, :);
    end
    
    plot(z/r, result.B_total/result.B_max, '-', ...
        'Color', line_color, 'LineWidth', 1.5, ...
        'DisplayName', 'Combined Field');
    
    % Formatting
    xlabel('Reduced distance from origin ($z/r$)', 'Interpreter', 'latex');
    ylabel('$B_z(0, z)/B_{z,max}(0, 0)$', 'Interpreter', 'latex');
    set(gca, 'TickLabelInterpreter', 'latex');
    
    if show_decomposition && size(result.B_individual, 1) > 1
        legend('Interpreter', 'latex', 'Location', 'Best');
    end
    
    % Create plots directory if it doesn't exist
    if ~exist('plots', 'dir')
        mkdir('plots');
    end
    
    % Save in both formats
    eps_filename = [plot_filename '.eps'];
    png_filename = [plot_filename '.png'];
    
    % Save EPS with high resolution
    print(fig, eps_filename, '-depsc', '-r300');
    
    % Save PNG with high resolution
    print(fig, png_filename, '-dpng', '-r300');
    
    % Close the figure immediately after saving
    close(fig);
    
    hold off;
end

%% Enhanced Plotting Function for Comparison Plots
function create_comparison_plot(z, r, results, plot_filename, plot_title, field_names, show_decomposition)
    if nargin < 7
        show_decomposition = false;
    end
    
    % Create figure with invisible property initially
    fig = figure('Units', 'centimeters', 'Position', [0, 0, 8, 6], 'Visible', 'off');
    hold on;
    
    % Special handling for plot 1f with comprehensive legend
    if show_decomposition && contains(plot_filename, 'plot1f')
        % Plot in specific order for a comprehensive legend
        for i = 1:length(field_names)
            config_name = field_names{i};
            result = results.(config_name);
            config = result.config;
            
            % Plot the combined field first
            plot(z/r, result.B_total/result.B_max, '-', ...
                'Color', config.color, 'LineWidth', 1, ...
                'DisplayName', config.name);
            
            % Plot decomposition for multi-coil systems with individual legend entries
            if size(result.B_individual, 1) > 1
                for j = 1:size(result.B_individual, 1)
                    % Create specific legend names for decomposition
                    decomp_name = sprintf('%s - Coil %d', config.name, j);
                    
                    plot(z/r, result.B_individual(j,:)/result.B_max, '--', ...
                        'Color', config.color, 'LineWidth', 0.5, ...
                        'DisplayName', decomp_name);
                end
            end
        end
    else
        % Standard plotting for other comparison plots
        % Counter for decomposition legend entries
        decomp_legend_added = false;
        
        for i = 1:length(field_names)
            config_name = field_names{i};
            result = results.(config_name);
            config = result.config;
            
            % Plot combined field
            plot(z/r, result.B_total/result.B_max, '-', ...
                'Color', config.color, 'LineWidth', 1, ...
                'DisplayName', config.name);
            
            % Plot decomposition for multi-coil systems if requested
            if show_decomposition && size(result.B_individual, 1) > 1
                for j = 1:size(result.B_individual, 1)
                    line_style = '--';
                    line_width = 0.5;
                    
                    % Add legend entry for decomposition lines
                    if ~decomp_legend_added
                        handle_vis = 'on';
                        display_name = 'Individual Coils';
                        decomp_legend_added = true;
                    else
                        handle_vis = 'off';
                        display_name = '';
                    end
                    
                    plot(z/r, result.B_individual(j,:)/result.B_max, line_style, ...
                        'Color', config.color, 'LineWidth', line_width, ...
                        'HandleVisibility', handle_vis, 'DisplayName', display_name);
                end
            end
        end
    end
    
    % Set labels based on plot type
    if show_decomposition
        xlabel('Reduced distance from origin ($z/r_n$)', 'Interpreter', 'latex');
        ylabel('$B_z(0, z)/B_z(0, 0)$', 'Interpreter', 'latex');
    else
        xlabel('Reduced distance from origin ($z/r_{sys}$)', 'Interpreter', 'latex');
        ylabel('$B_{z,sys}(0, z)/B_{z,sys}(0, 0)$', 'Interpreter', 'latex');
    end
    
    set(gca, 'TickLabelInterpreter', 'latex');
    legend('Interpreter', 'latex', 'Location', 'Best');
    
    % Create plots directory if it doesn't exist
    if ~exist('plots', 'dir')
        mkdir('plots');
    end
    
    % Save in both formats
    eps_filename = [plot_filename '.eps'];
    png_filename = [plot_filename '.png'];
    
    % Save EPS with high resolution
    print(fig, eps_filename, '-depsc', '-r300');
    
    % Save PNG with high resolution
    print(fig, png_filename, '-dpng', '-r300');
    
    % Close the figure immediately after saving
    close(fig);
    
    hold off;
end

%% Generate Individual Plots (Sequential Generation)
fprintf('\n=== Generating Individual Configuration Plots ===\n');

fprintf('Generating Single Coil plot...\n');
create_plot(z, r, results, 'single', 'plots/plot1a_v2_AWB');

fprintf('Generating Helmholtz Coil plot...\n');
create_plot(z, r, results, 'helmholtz', 'plots/plot1b_v2_AWB', true);

fprintf('Generating Maxwell Coil plot...\n');
create_plot(z, r, results, 'maxwell', 'plots/plot1c_v2_AWB', true);

fprintf('Generating Double-Pair Helmholtz Coil plot...\n');
create_plot(z, r, results, 'double_helmholtz', 'plots/plot1d_v2_AWB', true);

%% Generate Comparison Plots
fprintf('\n=== Generating Comparison Plots ===\n');

field_names = {'single', 'helmholtz', 'maxwell', 'double_helmholtz'};

fprintf('Generating Configuration Comparison plot...\n');
create_comparison_plot(z, r, results, 'plots/plot1e_v2_AWB', ...
    'All Configurations Comparison', field_names, false);

fprintf('Generating Detailed Comparison with Decomposition plot...\n');
create_comparison_plot(z, r, results, 'plots/plot1f_v2_AWB', ...
    'Detailed Comparison with Decomposition', field_names, true);

%% Performance Analysis
fprintf('\n=== Performance Analysis ===\n');
fprintf('Configuration\t\tMax Field (T)\tUniformity at z=0\n');
fprintf('%-20s\t%.6f\t%.6f\n', 'Single Coil', results.single.B_max, 1.0);
fprintf('%-20s\t%.6f\t%.6f\n', 'Helmholtz', results.helmholtz.B_max, 1.0);
fprintf('%-20s\t%.6f\t%.6f\n', 'Maxwell', results.maxwell.B_max, 1.0);
fprintf('%-20s\t%.6f\t%.6f\n', 'Double Helmholtz', results.double_helmholtz.B_max, 1.0);

%% Calculate uniformity metric (standard deviation in central region)
central_region = abs(z) <= 0.5*r;  % Central 50% of the range
fprintf('\n=== Field Uniformity Analysis ===\n');
for i = 1:length(field_names)
    config_name = field_names{i};
    result = results.(config_name);
    central_field = result.B_total(central_region);
    uniformity = std(central_field) / mean(central_field) * 100;  % Coefficient of variation
    fprintf('%-26s uniformity: %.2f%% variation\n', result.config.name, uniformity);
end

%% Ensure all figures are closed
close all;

%% Final Summary
fprintf('\n=== Magnetic Field Analyser v2.0 by Adil Wahab Bhatti ===\n');
fprintf('Advanced Analysis Complete! All plots saved to plots/ directory.\n');
fprintf('Generated plots in both EPS and PNG formats:\n');
fprintf('  • plot1a_v2_AWB (Single Coil)\n');
fprintf('  • plot1b_v2_AWB (Helmholtz Coil with decomposition)\n');
fprintf('  • plot1c_v2_AWB (Maxwell Coil with decomposition)\n');
fprintf('  • plot1d_v2_AWB (Double-Pair Helmholtz with decomposition)\n');
fprintf('  • plot1e_v2_AWB (All Configurations Comparison)\n');
fprintf('  • plot1f_v2_AWB (Detailed Comparison with Decomposition)\n');
fprintf('\nAll figure windows have been closed.\n');
