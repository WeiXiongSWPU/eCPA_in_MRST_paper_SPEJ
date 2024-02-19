classdef ECPABlackOilCapillaryPressure < StateFunction & SaturationProperty
    % Implementation of black-oil type capillary pressure
    properties
    end

    properties (Access = protected)
        surfaceTensionOW
        surfaceTensionOG
        porosityExponent
        permeabilityExponent
        permeabilityDirection
        pressureUnit = 1;
    end
    
    methods
        function prop = ECPABlackOilCapillaryPressure(model, varargin)
            prop = prop@StateFunction(model, varargin{:});
            prop = prop.dependsOn({'s', 'T'}, 'state');
            prop.label = 'p_{c}';
        end
        
        function pc = evaluateOnDomain(prop, model, state)
            [act, phInd] = model.getActivePhases();
            nph = sum(act);
            pc = cell(1, nph);
            
            f = model.fluid;
            JfuncActiveOW = prop.hasJFunctionScaler('OW');
            JfuncActiveOG = prop.hasJFunctionScaler('OG');
            if JfuncActiveOG || JfuncActiveOW
                ratio = prop.getJFunctionStaticRatio(model);
            end

            if model.water && model.oil && isfield(f, 'pcOW')
                sW = model.getProp(state, 'sw');
                if prop.scalingActive
                    pts = model.rock.krscale.drainage;
                    reg = prop.regions;
                    [get, ~, U, L] = SaturationProperty.getSatPointPicker(f, pts, reg, prop.cell_subset);
                    [swcon, SWCON] = get('w', L);
                    [swmax, SWMAX] = get('w', U);
                    sW = swcon + (sW - SWCON).*(swmax - swcon)./(SWMAX - SWCON);
                end
                pcow = prop.evaluateFluid(model, 'pcOW', sW);
                if isfield(state, 'pcowScale')
                    pcow = pcow.*state.pcowScale;
                    assert(~JfuncActiveOW, 'Cannot both have initial water pc scale and Jfunc scaling.');
                elseif JfuncActiveOW
                    sow = prop.surfaceTensionOW;
                    pcow = sow.*ratio.*pcow;
                end
                % Note sign! Water is always first
                pc{phInd == 1} = -pcow;
            end
            
            if model.gas && model.oil && isfield(f, 'pcOG')
                % Oil gas capillary pressure
                sG = model.getProp(state, 'sg');
                pcog = prop.evaluateFluid(model, 'pcOG', sG);
                if JfuncActiveOG
                    sog = prop.surfaceTensionOG;
                    pcog = sog.*ratio.*pcog;
                end
                pc{phInd == 3} = pcog;
            end
            
            if ~model.oil && isfield(f, 'pcWG')
                % Water-gas capillary pressure
                sW = model.getProp(state, 'sw');
                if isfield(model.rock, 'BC')
                    Pce = model.rock.BC.Pce;
                    sWi = model.rock.BC.sWi;
                    pc{2} = Pce .* ((sW-sWi) ./ (1-sWi)).^(-0.5);
                elseif isfield(model.rock, 'VG')
                    Pce = model.rock.VG.Pce;
                    sWi = model.rock.VG.sWi;
                    snt = model.rock.VG.snt;
                    select = value(sW)>(1-snt);
                    pc{2} = Pce .* ((sW-sWi) ./ (1-sWi)).^(-0.5);
                    if any(select)
                        pc{2}(select) = Pce./snt .* ((1-snt-sWi)./(1-sWi)).^(-0.5).*(1-sW(select));
                    end
%                 elseif isfield(model.rock, 'Silin')
%                     A = model.rock.Silin.A;
%                     B = model.rock.Silin.B;
%                     alpha1 = -model.rock.Silin.alpha1;
%                     alpha2 = model.rock.Silin.alpha2;
%                     sWi = model.rock.Silin.sWi;
%                     sWs = (sW-sWi) ./ (1-sWi);
%                     pc{2} = A.*(sWs.^alpha1-1)+B.*(1-sWs.^alpha2).^(1./alpha2);
%                     pc = ensureMaxDerivatives(prop, pc);
%                 elseif isfield(model.rock, 'BCIFT')
%                     T = model.getProp(state, 'T');
%                     L = model.getProp(state, 'L');
%                     x = state.x;
%                     IFT_pure = eCPApureInterfacialTension(model.EOSModel, T);
%                     IFT = eCPAInterfacialTension(model.EOSModel, T, x, IFT_pure, L);
%                     
%                     phi = model.rock.poro;
%                     k = model.rock.perm;
%                     ratio = phi./k;
%                     
%                     Pce = model.rock.BCIFT.Pce;
%                     sWi = model.rock.BCIFT.sWi;
%                     pc{2} = Pce .*((sW-sWi) ./ (1-sWi)).^(-0.5) .* IFT .* ratio.^0.5;
                elseif isfield(model.rock, 'VGIFT')
                    T = model.getProp(state, 'T');
                    L = model.getProp(state, 'L');
                    x = state.x;
                    IFT_pure = eCPApureInterfacialTension(model.EOSModel, T);
                    IFT = eCPAInterfacialTension(model.EOSModel, T, x, IFT_pure, L);
                    
                    phi = model.rock.poro;
                    k = model.rock.perm;
                    ratio = phi./k(:,1);
                    
                    Pce = model.rock.VGIFT.Pce;
                    sWi = model.rock.VGIFT.sWi;
                    snt = model.rock.VGIFT.snt;
                    pc{2} = Pce .*((sW-sWi) ./ (1-sWi)).^(-0.5) .* IFT .* ratio.^0.5;
                    
                    select = value(sW)>(1-snt);
                    if any(select)
                        pc{2}(select) = Pce./snt .* ((1-snt-sWi)./(1-sWi)).^(-0.5).*(1-sW(select)).* IFT(select) .* ratio(select).^0.5;
                    end
                else
                    sG = model.getProp(state, 'sg');
                    pc{phInd == 3} = prop.evaluateFluid(model, 'pcWG', sG);
                end
            end
        end
        
%         function Pc = ensureMaxDerivatives(prop, Pc)
%             nc = numel(Pc);
%             for i = 1:nc
%                 if ~isa(Pc{i}, 'ADI') && ~isempty(Pc{i})
%                     return
%                 end
%             end
%             der = 1e3;
%             if numel(der) == 1
%                 der = repmat(der, 1, nc);
%             end
%             isDiag = isa(prop.AutoDiffBackend, 'DiagonalAutoDiffBackend');
%             if isDiag
%                 rowMajor = prop.AutoDiffBackend.rowMajor;
%                 for c = 1:nc
%                     m = Pc{c};
%                     if isnumeric(m) || size(m.jac{1}.diagonal, 2 - rowMajor) < c
%                         continue
%                     end
%                     d = der(c);
%                     if rowMajor
%                         bad = abs(m.jac{1}.diagonal(c, :)) > d;
%                         m.jac{1}.diagonal(c, bad) = d;
%                     else
%                         bad = abs(m.jac{1}.diagonal(:, c)) > d;
%                         m.jac{1}.diagonal(bad, c) = d;
%                     end
%                     Pc{c} = m;
%                 end
%             else
%                 for c = 1:nc
%                     m = Pc{c};
%                     d = der(c);
%                     if isnumeric(m)
%                         continue;
%                     end
%                     if numel(m.jac) < c
%                         % Not matching number of derivatives - stop early
%                         break;
%                     end
%                     J = m.jac{c};
%                     [n, l] = size(J);
%                     if n ~= l
%                         continue;
%                     end
%                     diagonal = diag(J);
%                     bad = abs(diagonal) > d;
%                     if any(bad)
%                         m.jac{c} = m.jac{c} + sparse(1:n, 1:n, bad.*d, n, n);
%                     end
%                     Pc{c} = m;
%                 end
%             end
%         end
        
        function anyPresent = pcPresent(prop, model)
            f = model.fluid;
            anyPresent = isfield(f, 'pcOW') || isfield(f, 'pcOG') || isfield(f, 'pcWG');
        end
        
        function property = subset(property, subs)
            property = subset@StateFunction(property, subs);
            property.cell_subset = subs;
        end

        function prop = setJFunctionConstants(prop, poroexp, permexp, permdir, pressureunit)
            if nargin > 4
                prop.pressureUnit = pressureunit;
            end
            prop.porosityExponent = poroexp;
            prop.permeabilityExponent = permexp;
            
            if ischar(permdir)
                newdir = zeros(1, numel(permdir));
                for i = 1:numel(permdir)
                    newdir(i) = find('xyz' == lower(permdir(i)));
                end
                permdir = newdir;
            end
            prop.permeabilityDirection = permdir; %
        end
        
        function prop = setSurfaceTension(prop, value, fluidpair)
            switch lower(fluidpair)
                case 'ow'
                    prop.surfaceTensionOW = value;
                case 'og'
                    prop.surfaceTensionOG = value;
                otherwise
                    error('Unsupported pair %s', fluidpair);
            end
        end
        
        function present = hasJFunctionScaler(prop, phasepair)
            present = prop.scalingActive && ~isempty(prop.getSurfaceTension(phasepair));
        end
        
        function st = getSurfaceTension(prop, phasepair)
            nm = ['surfaceTension', upper(phasepair)];
            st = prop.(nm);
        end
        
        
        function ratio = getJFunctionStaticRatio(prop, model)
            phi = model.rock.poro(prop.cell_subset);
            k = model.rock.perm(prop.cell_subset, prop.permeabilityDirection);
            k = sum(k, 2)./size(k, 2);
            % Apply exponents
            k = k.^prop.permeabilityExponent;
            phi = phi.^prop.porosityExponent;
            ratio = phi./k;
            ratio = ratio./prop.pressureUnit; % Account for unit for table
        end
        
    end
end

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
