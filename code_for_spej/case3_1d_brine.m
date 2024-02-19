%% Simulate an one-dimensional example for brine-CO2
% This script contains the example case 3 for the paper
% "Compositional simulation for carbon storage in porous media using
%  an electrolyte association equation of state"
%
% To demonstrates the difference in practice between the two compositional
% formulations (Natural Variable and Overall Composition Formulations)
% implemented in MRST, and also show to set up the corresponding solvers,
% we consider a simplified two-component water-CO2 system
% posed inside 100 grid cells.
mrstModule add compositional ad-core ad-props
mrstVerbose on

%% Construct the two simulation models
% Grid and petrophysical data in 100 cells
G    = computeGeometry(cartGrid(100));      % 100 cells 1 m^3 grid
rock = makeRock(G, 0.5*darcy, 0.5);

% Start with fluid properties from a standard black-oil, water-gas model
% and add extra mixture properties to turn this into a brine-CO2 model.
f = initSimpleADIFluid('phases', 'wg', 'blackoil', false, 'rho', [1000, 700]);
mixture = ECPATableCompositionalMixture({'Water','Carbondioxide', 'Na+','Cl-'});

% Construct models for both formulations. Same input arguments
arg = {G, rock, f, ...                              % Standard arguments
       mixture,...                                  % Compositional mixture
       'water', true, 'oil', false, 'gas', true,... % Water-Gas system
       'liquidPhase', 'W', 'vaporPhase', 'G'};      % Water=liquid, gas=vapor
diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'rowMajor', true);
overall = ECPAGenericOverallCompositionModel(arg{:}, 'AutoDiffBackend', diagonal_backend);   % Overall mole fractions
natural = ECPAGenericNaturalVariablesModel(arg{:}, 'AutoDiffBackend', diagonal_backend);     % Natural variables

% Validate both models to initialize the necessary state function groups
overall = overall.validateModel();
natural = natural.validateModel();

%% Set up boundary conditions and initial conditions
% We start by the boundary conditions. Here, we use standard routines from
% MRST, but must remember to also impose component specification at all
% boundaries
p = 50*barsa; T = 273.15 + 30; s = []; z = [0.917,0.05,0.0165,0.0165];     % p, T, z, s
bc = fluxside([], G, 'xmin', 1/day, 'sat', [0, 1]);    % Flux
bc = pside(bc, G, 'xmax', p, 'sat', [0, 1]);           % Standard bc
bc.components = repmat([0, 1,0,0], numel(bc.face), 1);     % Boundary z

% We first set the initial state from the prescribed parameters, validate
% the model (which will add data structures for any wells, etc), and then
% perform a flash to ensure that the corresponding initial guess exists at
% equilibrium conditions
state0 = eCPAinitCompositionalState(overall, p, T, s, z);
state0 = overall.validateState(state0);
initGuess = state0;
initGuess = overall.computeFlash(initGuess);

%% Adjust solver settings and construct nonlinear solver
% To improve the nonlinear solver we impose limits on changes during the
% nonlinear iterations. With natural variables, we limit phase mole
% fractions and saturation, whereas for the overall composition we limit
% the overall mole fraction
    % change
if ~exist('maxChange', 'var')
    maxChange = 0.1;
end
if isfinite(maxChange)
    overall.dzMaxAbs = maxChange;
    natural.dzMaxAbs = maxChange;
    natural.dsMaxAbs = maxChange;
end

% Set up nonlinear solver with high report level to output intermediate
% states
nls = NonLinearSolver('reportLevel', 3, 'MaxIterations', 1000);


%% Compute a single time step
dt = 100*day;
[solNat, reportsNat] = nls.solveTimestep(state0, dt, natural, ...
                            'bc', bc, 'initialGuess', initGuess);
[solMole, reportsMole] = nls.solveTimestep(state0, dt, overall, ...
                            'bc', bc, 'initialGuess', initGuess);

%% Extract data from the output states
getStates = @(reports) cellfun(@(x) x.state, reports.StepReports{1}.NonlinearReport, 'UniformOutput', false);
% Get the outputs
getZ = @(states) [initGuess.components(1); cellfun(@(x) x.components(1), states)];
getS = @(states) [initGuess.s(1); cellfun(@(x) x.s(1), states)];
getP = @(states) [initGuess.pressure(1); cellfun(@(x) x.pressure(1), states)] - p;
natStates = getStates(reportsNat);
moleStates = getStates(reportsMole);

zn = getZ(natStates);
sn = getS(natStates);
dpn = getP(natStates);

zm = getZ(moleStates);
sm = getS(moleStates);
dpm = getP(moleStates);
%% Plot the convergence histories for both models
% We make two different plots of the convergence: in (z_w,dp) space and in
% (S_l,dp) space. These plots are imposed on top of a contour map of S.
plotConvergence(overall, p, T, zn, sn, dpn, zm, sm, dpm, maxChange)
plotSolution(solNat, solMole)

%% Plotting functions
function plotSolution(solNat, solMole)
    nc = numel(solNat.T);
    figure
    h1 = plot(1:nc, solNat.pressure, '-o', 'color', 'k', 'linewidth', 1.2, ...
        'MarkerFaceColor', 'r', 'MarkerSize', 8); hold on
    h2 = plot(1:nc, solMole.pressure, '-o', 'color', 'b', 'linewidth', 1.2, ...
        'MarkerFaceColor', 'y');
    legend([h1; h2], 'Natural variables', 'Overall composition', ...
         'location', 'southeast')
    ylabel('\Delta p [Pa]')
    xlabel('m')
    figure
    h1 = plot(1:nc, solNat.s(:,1), '-o', 'color', 'k', 'linewidth', 1.2, ...
        'MarkerFaceColor', 'r', 'MarkerSize', 8); hold on
    h2 = plot(1:nc, solMole.s(:,1), '-o', 'color', 'b', 'linewidth', 1.2, ...
        'MarkerFaceColor', 'y');
    legend([h1; h2], 'Natural variables', 'Overall composition', ...
         'location', 'southeast')
    ylabel('s')
    xlabel('m')
    figure
    h1 = plot(1:nc, solNat.components(:,1), '-o', 'color', 'k', 'linewidth', 1.2, ...
        'MarkerFaceColor', 'r', 'MarkerSize', 8); hold on
    h2 = plot(1:nc, solMole.components(:,1), '-o', 'color', 'b', 'linewidth', 1.2, ...
        'MarkerFaceColor', 'y');
    legend([h1; h2], 'Natural variables', 'Overall composition', ...
         'location', 'southeast')
    ylabel('z')
    xlabel('m')
end

function plotConvergence(overall, p, T, zn, sn, dpn, zm, sm, dpm, maxChange)
    ns = 50;
    ps = linspace(min(min(dpn), min(dpm)), 1.1*max(max(dpn), max(dpm)), ns);
    zs = linspace(0.01, 1, ns);

    [DP, Z] = meshgrid(ps, zs);
    zz = reshape(Z, [], 1);

    [L, ~, ~, Z_L, Z_V] = standaloneFlash(DP(:) + p, T, [zz.*55.51./(55.51+2), 1 - zz,zz./(55.51+2),zz./(55.51+2)], overall.EOSModel);
    S  = reshape(L.*Z_L./(L.*Z_L + (1-L).*Z_V), ns, ns);

    for i = 1:2
        if i == 1
            xn = zn; xm = zm; XX = Z; l = 'z_{water}';
        else
            xn = sn; xm = sm; XX = S; l = 'S_L';
        end
        figure; clf; hold on
        contourf(XX, DP, S, 10);
        c = lines(4);
        c1 = c(2, :);
        c2 = c(3, :);
        if maxChange < 0.1
            style = '-';
        else
            style = '-o';
        end
        h1 = plot(xn, dpn, style, 'color', c1, 'linewidth', 1.2, ...
            'MarkerFaceColor', 'r', 'MarkerSize', 8); hold on
        h2 = plot(xm, dpm, style, 'color', c2, 'linewidth', 1.2, ...
            'MarkerFaceColor', 'y');
        h3 = plot(xm(end), dpm(end), 'x', 'color', [1, 1, 1]*0.7, ...
            'MarkerSize', 10, 'linewidth', 1.2);
        legend([h1; h2; h3], 'Natural variables', 'Overall composition', ...
            'Solution', 'location', 'southeast')
        ylabel('\Delta p [Pa]')
        xlabel(l)

        cmap = interp1([1; 100], [.1, .1, 1; 0.9, 0.9, 1], linspace(1, 100, 100));
        colormap(flipud(cmap))
    end
end

%%
%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
