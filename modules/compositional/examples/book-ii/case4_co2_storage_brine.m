%% Simulation CO2 storage in saline aquifers
% Geological storage of CO2 into saline aquifers consists of two periods
% which include injection and post injection periods. The model developed
% is used to study CO2 dissolution into formation brines during the
% injection period and using an idealized 2D geometry for the aquifer.
mrstModule add ad-core ad-props mrst-gui compositional deckformat
gravity reset on

% The numerical model consists of 99 grid blocks in horizontal direction
% and 39 grid blocks in vertical direction.The grid blocks in horizontal
% direction are refined around the injection well, with a minimum of one
% meter close to the injection well and a maximum of 99 meters at the outer
% boundary, which is located 5 km away from the injection well. The grid
% blocks in the vertical direction have the same size of one meter.
dx = 0:99;
nx = cumsum(dx);
ny = 0:1;
nz = 0:39;

G = tensorGrid(nx, ny, nz);
G = computeGeometry(G);
plotGrid(G); view(3); axis tight

% Assume an isotropic and homogenous aquifer with a permeability of 100 mD,
% porosity of 0.18
rock = makeRock(G, 100*milli*darcy, 0.18);

% Rock compressibility is 8e-5/bar
f = initSimpleADIFluid('phases', 'wg', 'blackoil', false, 'rho', [1000, 700], 'cR', 8e-5/barsa, 'n', [4, 2]);
ECPAmixture = ECPATableCompositionalMixture({'Water','Carbondioxide','Na+','Cl-'});

% Construct models for both formulations. Same input arguments
ECPAarg = {G, rock, f, ...                              % Standard arguments
       ECPAmixture,...                                  % Compositional mixture
       'water', true, 'oil', false, 'gas', true,...     % Water-Gas system
       'liquidPhase', 'W', 'vaporPhase', 'G'};          % Water=liquid, gas=vapor
diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'rowMajor', true);
ECPAoverall = ECPAGenericOverallCompositionModel(ECPAarg{:}, 'AutoDiffBackend', diagonal_backend);   % Overall mole fractions
ECPAnatural = ECPAGenericNaturalVariablesModel(ECPAarg{:}, 'AutoDiffBackend', diagonal_backend);     % Natural variables
ECPAoverall = imposeRelpermScaling(ECPAoverall, 'SWCR', 0.3, 'SGU', 0.7);
ECPAnatural = imposeRelpermScaling(ECPAnatural, 'SWCR', 0.3, 'SGU', 0.7);
% Validate both models to initialize the necessary state function groups
ECPAoverall = ECPAoverall.validateModel();
ECPAnatural = ECPAnatural.validateModel();

%% Set up BC + initial conditions/initial guess
% The impermeable top layer of the aquifer is located at a depth of 1200 m
% with a corresponding pressure and temperature of 12 MPa and 38 ¡æ. brine
% have 40000 ppm salinity.
p = 120*barsa; T = 38+273.15; z = [0.976,0,0.012,0.012];     % p, T, z

% The simulations are continued for up to 150 years after CO2 injection has
% stopped.
totTime = 150*year;

% A constant pressure boundary condition is imposed on the far boundary at
% the right hand side of the domain, 5 km away from the injector, to
% simulate the aquifer outflow and avoid over pressurization of the
% aquifer.
bc = pside([], G, 'xmax', p, 'sat', [0, 1]);       % Standard bc
bc.components = repmat(z, numel(bc.face), 1); % Boundary z

% CO2 is injected at a rate of 1000 m3/day per meter width of the aquifer
% for 2 years into a well located at the left hand side of the domain which
% is completed in the bottom grid block.
W = addWell([], G, rock, 3763, 'name', 'inj', 'type', 'rate', 'Compi', ...
[0 1], 'val', 1000/day, 'components', [0 1 0 0]);

% Plot well
show = true([G.cells.num, 1]);
cellInx = 1:99:3861;
show(cellInx) = false;
g = norm(gravity);
[z_0, z_max] = deal(0, max(G.cells.centroids(:,3)));
equil  = ode23(@(z,p) g .* 1000, [z_0, z_max], p);
p_init = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);  clear equil
plotCellData(G, convertTo(p_init, barsa), show, 'EdgeColor', 'k')
plotWell(G, W, 'height', 10)
view(-125, 20), camproj perspective

% Initialize the problem
s0 = [1, 0];
state0 = initResSol(G, p, s0);
state0.T = repmat(T, G.cells.num, 1);
state0.components = repmat(z, G.cells.num, 1);

% CO2 is injected at a rate of 1000 m3/day per meter width of the aquifer
% for 2 years into a well located at the left hand side of the domain which
% is completed in the bottom grid block.
dt1 = rampupTimesteps(2*year, 75*day, 12);
dt2 = rampupTimesteps(totTime-2*year, 90*day, 12);
dt=[dt1;dt2];

schedule = simpleSchedule(dt,'W', W, 'bc', bc);
schedule.step.control(numel(dt1)+1:end) = 2;
schedule.control(2) = schedule.control(1);
schedule.control(2).W.val = 0;

% Set up nonlinear solver with high report level to output intermediate
% states
nls = NonLinearSolver('useRelaxation', true);
% [~, ECPAstatesn, ECPAreportn] = simulateScheduleAD(state0, ECPAnatural, schedule, 'nonlinearsolver', nls);
[~, ECPAstateso, ECPAreporto] = simulateScheduleAD(state0, ECPAoverall, schedule, 'nonlinearsolver', nls);

%% Launch interactive plotting
figure;
plotToolbar(G, ECPAstateso)
view(0,0)
colorbar

f=zeros(1,numel(ECPAstateso));
for i = 1:numel(ECPAstateso)
    L =ECPAstateso{i}.L;
    y =ECPAstateso{i}.y(:,1);
    n=y.*(1-L);
    X = sum(n);
    f(i) = X;
end
figure
plot(cumsum(dt./year), f)


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
