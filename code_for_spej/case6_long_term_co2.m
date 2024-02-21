%% Benchmark for long-term CO2 sequestration simulations
mrstModule add ad-core ad-props mrst-gui compositional deckformat linearsolvers
gravity reset on

nx = 0:403.86:16154.4;
ny = 0:403.86:16154.4;
nz = 0:7.62:304.8;
dz = (161.544:-4.0386:0);
dz = meshgrid(dz ,dz);

G = tensorGrid(nx, ny, nz,'depthz', dz');
G = computeGeometry(G);
plotGrid(G); view(3); axis tight

rock = makeRock(G, [100,100,100]*milli*darcy, 0.25);
[nx, ny, nz] = deal(G.cartDims(1), G.cartDims(2), G.cartDims(3));
nf = nx*ny;
vhr = 0.001;
rock.perm(1:4*nf,:) = repmat([89,89,89*vhr]*milli*darcy,4*nf,1);
rock.perm(1+4*nf:8*nf,:) = repmat([65,65,65*vhr]*milli*darcy,4*nf,1);
rock.perm(1+8*nf:12*nf,:) = repmat([46,46,46*vhr]*milli*darcy,4*nf,1);
rock.perm(1+12*nf:16*nf,:) = repmat([30,30,30*vhr]*milli*darcy,4*nf,1);
rock.perm(1+16*nf:20*nf,:) = repmat([15,15,15*vhr]*milli*darcy,4*nf,1);
rock.perm(1+20*nf:24*nf,:) = repmat([120,120,120*vhr]*milli*darcy,4*nf,1);
rock.perm(1+24*nf:28*nf,:) = repmat([165,165,165*vhr]*milli*darcy,4*nf,1);
rock.perm(1+28*nf:32*nf,:) = repmat([235,235,235*vhr]*milli*darcy,4*nf,1);
rock.perm(1+32*nf:36*nf,:) = repmat([840,840,840*vhr]*milli*darcy,4*nf,1);
rock.perm(1+36*nf:40*nf,:) = repmat([370,370,370*vhr]*milli*darcy,4*nf,1);

f = initSimpleADIFluid('phases', 'wg', 'blackoil', false, 'rho', [1000, 700], 'n', [2, 2]);
ECPAmixture = ECPATableCompositionalMixture({'Water','Carbondioxide'});

% Construct models for both formulations. Same input arguments
ECPAarg = {G, rock, f, ...                              % Standard arguments
       ECPAmixture,...                                  % Compositional mixture
       'water', true, 'oil', false, 'gas', true,...     % Water-Gas system
       'liquidPhase', 'W', 'vaporPhase', 'G'};          % Water=liquid, gas=vapor
diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'rowMajor', true);
ECPAoverall = ECPAGenericOverallCompositionModel(ECPAarg{:}, 'AutoDiffBackend', diagonal_backend);   % Overall mole fractions
ECPAnatural = ECPAGenericNaturalVariablesModel(ECPAarg{:}, 'AutoDiffBackend', diagonal_backend);     % Natural variables
ECPAoverall = imposeRelpermScaling(ECPAoverall, 'SWCR', 0.25,  'KRW', 0.334, 'SGCR', 0.25,'KRG',1,'SGU',0.75,'SWU',0.75);
ECPAnatural = imposeRelpermScaling(ECPAnatural, 'SWCR', 0.25,  'KRW', 0.334, 'SGCR', 0.25,'KRG',1,'SGU',0.75,'SWU',0.75);
% Validate both models to initialize the necessary state function groups
ECPAoverall = ECPAoverall.validateModel();
ECPAnatural = ECPAnatural.validateModel();

%% Set up BC + initial conditions/initial guess
p0 = 156*barsa; T = 60+273.15; z = [1,0];     % p, T, z

eCPA = ECPAEquationOfStateModel([], ECPAmixture, 'eCPA');
[~, ~, ~,~, ~,rho0] = eCPAstandaloneFlash(p0, T, z, eCPA);  
g = norm(gravity);
[z_0, z_max] = deal(0, max(G.cells.centroids(:,3)));
equil  = ode23(@(z,p) g .* rho0, [z_0, z_max], p0);
p = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);  clear equil

totTime = 1000*year;

% constant pressure boundary conditions
bc = pside([], G, 'xmax', p(G.cells.centroids(:,1)==15952.47), 'sat', [1, 0]);       % Standard bc
bc = pside(bc, G, 'xmin', p(G.cells.centroids(:,1)<220), 'sat', [1, 0]);       % Standard bc
bc = pside(bc, G, 'ymax', p(G.cells.centroids(:,2)==15952.47), 'sat', [1, 0]);       % Standard bc
bc = pside(bc, G, 'ymin', p(G.cells.centroids(:,2)<220), 'sat', [1, 0]);       % Standard bc
bc.components = repmat(z, numel(bc.face), 1);      % Boundary z

% CO2 is injected at a rate of 960000 t/year for 50 years
[~, ~, ~,~, ~,rho] = eCPAstandaloneFlash(1*barsa, 298.15, [0, 1], eCPA);  
W = addWell([], G, rock, 63170, 'name', 'inj', 'type', 'rate', 'Compi', ...
[0 1], 'val', 960000*1000/rho/year, 'components', [0, 1]);

% Plot well
show = true([G.cells.num, 1]);
cellInx = 770:1600:63170;
show(cellInx) = false;
plotCellData(G, convertTo(p, barsa), show, 'EdgeColor', 'k')
plotWell(G, W, 'height', 10)
view(0, 0), camproj perspective

% Initialize the problem
s0 = [1, 0];
state0 = initResSol(G, p, s0);
state0.T = repmat(T, G.cells.num, 1);
state0.components = repmat(z, G.cells.num, 1);

dt1 = rampupTimesteps(50*year, 1*year, 10);
dt2 = rampupTimesteps(totTime-50*year, 1*year, 10);
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

figure
plotToolbar(G, ECPAstateso)
view(0,0)
colorbar


%%
%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
