function [IFT, rho_L, rho_V] = eCPApureInterfacialTension(eos, T)
% Estimate interfacial tension of pure substances

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

rhoc = 1./eos.ECPACompositionalMixture.Vcrit;
names = eos.ECPACompositionalMixture.names;
nmole = eos.getNumberOfMolecules();
parachor = eos.ECPACompositionalMixture.parachor;
Tc = eos.ECPACompositionalMixture.Tcrit;
Tr = T ./ Tc;

[a, b, m, n, p] = getIFTparameters({names{1:nmole}});
tao = 1 - Tr(:,1:nmole);
beta = 1 + tao.^b;
rho_L = rhoc(1:nmole).*exp(a.*(beta-exp(1-beta)));

alpha = exp(tao.^(1/3) +tao.^0.5 +tao +tao.^m);
rho_V = rhoc(1:nmole).*exp(p.*(alpha.^n-exp(1-alpha)));

IFT=1e-3 * (parachor(1:nmole).*(rho_L-rho_V ).*1e-6).^4;

bad = T > Tc;
IFT(bad) = 0;
rho_L(bad) = nan;
rho_V(bad) = nan;
end

function [a, b, m, n, p] = getIFTparameters(names)
fluids = IFTPropFluidsStructs();
validChoices = {fluids.name};
ok = ismember(lower(names), lower(validChoices));
if ~all(ok)
    s ='Unable to create fluid. The following names were not known to CoolProps table: ';
    msg = [s, sprintf('%s ', names{~ok})];
    error(msg);
end
ncomp = numel(names);
[a, b, m, n, p] = deal(zeros(1, ncomp));
for i = 1:ncomp
    isF = strcmpi(names{i}, validChoices);
    str = fluids(isF);
    a(i) = str.a;
    b(i) = str.b;
    m(i) = str.m;
    n(i) = str.n;
    p(i) = str.p;
end
end

function fluids = IFTPropFluidsStructs()
fluids = [...
    struct('name', 'Water', 'a', 0.828949, 'b', 0.287388, 'm', 2.3609558, 'n', 1.0916682, 'p', -0.7452828), ...
    struct('name', 'CarbonDioxide', 'a', 0.8003152, 'b', 0.3245757, 'm', 2.4686277, 'n', 1.1345838, 'p', -0.6240188), ...
    struct('name', 'Methane', 'a', 0.7482718, 'b', 0.3283942, 'm', 2.9407241, 'n', 1.1012273, 'p', -0.5697413), ...
    struct('name', 'HydrogenSulfide', 'a', 0.7874, 'b', 0.3406, 'm', 2.988, 'n', 1.15, 'p', -0.568), ...
    struct('name', 'Nitrogen', 'a', 0.7572631, 'b', 0.3312663, 'm', 2.8539512, 'n', 1.1131452, 'p', -0.5673304), ...
    struct('name', 'Oxygen', 'a', 0.7544534, 'b', 0.3311013, 'm', 2.9401387, 'n', 1.111504, 'p', -0.5649451), ...
    struct('name', 'Argon', 'a', 0.7489134, 'b', 0.3317864, 'm', 3.1601166, 'n', 1.1183779, 'p', -0.5515808), ...
    struct('name', 'Krypton', 'a', 0.7592822, 'b', 0.3288551, 'm', 3.4377347, 'n', 1.1459094, 'p', -0.5312377), ...
    struct('name', 'Ethane', 'a', 0.7674907, 'b', 0.3259824, 'm', 2.6379996, 'n', 1.1048142, 'p', -0.6018817), ...
    struct('name', 'n-Propane', 'a', 0.7828663, 'b', 0.3191917, 'm', 2.5491601, 'n', 1.1169908, 'p', -0.6110244), ...
    struct('name', 'n-Butane', 'a', 0.7858455, 'b', 0.318336, 'm', 2.470393, 'n', 1.1187258, 'p', -0.629651003), ...
    struct('name', 'n-Pentane', 'a', 0.7912326, 'b', 0.3154526, 'm', 2.4715148, 'n', 1.132658, 'p', -0.6412635), ...
    struct('name', 'n-Hexane', 'a', 0.7993375, 'b', 0.3119135, 'm', 2.5036259, 'n', 1.1549903, 'p', -0.6410813), ...
    struct('name', 'n-Heptane', 'a', 0.8148013, 'b', 0.311786, 'm', 2.547772, 'n', 1.1770021, 'p', -0.6419602), ...
    struct('name', 'Xenon', 'a', 0.7571569, 'b', 0.3203486, 'm', 3.2336283, 'n', 1.1272991, 'p', -0.5475391), ...
    struct('name', 'Neon', 'a', 0.7440979, 'b', 0.3421012, 'm', 4.4960428, 'n', 1.1861523, 'p', -0.485720905), ...
    struct('name', 'Fluorine', 'a', 0.7436589, 'b', 0.3499347, 'm', 2.5237557, 'n', 1.0739354, 'p', -0.6052029), ...
    struct('name', 'CarbonMonoxide', 'a', 0.7683354, 'b', 0.3320236, 'm', 2.1804625, 'n', 1.0460487, 'p', -0.6091577), ...
    struct('name', 'Ethylene', 'a', 0.771377, 'b', 0.3285082, 'm', 2.6032071, 'n', 1.1016984, 'p', -0.5971817), ...
    struct('name', 'Ammonia', 'a', 0.8598304, 'b', 0.309669, 'm', 2.9025748, 'n', 1.1747326, 'p', -0.6213074), ...
    struct('name', 'Propylene', 'a', 0.8019922, 'b', 0.3249273, 'm', 3.1010584, 'n', 1.187285, 'p', -0.5517088), ...
    struct('name', 'IsoButane', 'a', 0.789092, 'b', 0.3196161, 'm', 2.5860317, 'n', 1.1380881, 'p', -0.6069818), ...
    struct('name', 'R22', 'a', 0.7974639, 'b', 0.3234936, 'm', 2.447419, 'n', 1.1245385, 'p', -0.6346948), ...
    struct('name', 'R143a', 'a', 0.8144805, 'b', 0.3185632, 'm', 2.4390945, 'n', 1.1210065, 'p', -0.6625750), ...
    struct('name', 'R152A', 'a', 0.8210739, 'b', 0.3183938, 'm', 2.4036109, 'n', 1.1186849, 'p', -0.6709599), ...
    struct('name', 'R32', 'a', 0.8447267, 'b', 0.3189232, 'm', 2.4973193, 'n', 1.1207203, 'p', -0.6786065), ...
    struct('name', 'R123', 'a', 0.7993085, 'b', 0.3171035, 'm', 2.3683423, 'n', 1.137429, 'p', -0.6444662), ...
    struct('name', 'R124', 'a', 0.7973294, 'b', 0.3200154, 'm', 2.4354428, 'n', 1.1506277, 'p', -0.6362176), ...
    struct('name', 'R125', 'a', 0.7990249, 'b', 0.3167709, 'm', 2.3507937, 'n', 1.1406387, 'p', -0.6493031), ...
    struct('name', 'R134a', 'a', 0.8138468, 'b', 0.3182686, 'm', 2.420336, 'n', 1.1516298, 'p', -0.6546052), ...
    struct('name', 'NitrogenTrifluoride', 'a', 0.7924894, 'b', 0.3082062, 'm', 2.7962186, 'n', 1.1483585, 'p', -0.5712755)];

end