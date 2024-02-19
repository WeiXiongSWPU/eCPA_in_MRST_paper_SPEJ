classdef ECPAPhaseMixingCoefficientsLV < StateFunction
    % Mixing coefficients for liquid-vapor system
    properties
        useCompactEvaluation = true;
    end
    
    methods
        function gp = ECPAPhaseMixingCoefficientsLV(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'T'}, 'state');
            gp = gp.dependsOn('ComponentPhaseMoleFractions');
            gp.label = 'acti A B a bic Tr';
        end

        function v = evaluateOnDomain(prop, model, state)
            eos = model.EOSModel;
            T = model.getProps(state, 'temperature');
            xy = prop.getEvaluatedDependencies(state, 'ComponentPhaseMoleFractions');
            L_ix = model.getLiquidIndex();
            V_ix = model.getVaporIndex();
            isEoS = model.getEoSComponentMask();
            nph = size(xy, 2);
            v = cell(1, nph);
            
            x = xy(isEoS, L_ix);
            y = xy(isEoS, V_ix);
            
            [bic, a, Tr] = eos.getBic(T, iscell(x));
            [A_L, B_L, acti_L] = eos.getPhaseMixCoefficients(x, T, bic, a);
            twoPhase = model.getTwoPhaseFlag(state);
            if prop.useCompactEvaluation && iscell(y) && ~all(twoPhase)
                A_V = A_L;
                B_V = B_L;
                acti_V = acti_L;
                if any(twoPhase)
                    bic_2ph = cellfun(@(x) x(twoPhase), bic, 'UniformOutput', false);
                    a_2ph = cellfun(@(x) x(twoPhase), a, 'UniformOutput', false);
                    y_2ph = cellfun(@(x) x(twoPhase), y, 'UniformOutput', false);
                    [A_V(twoPhase), B_V(twoPhase), acti] = eos.getPhaseMixCoefficients(y_2ph, T(twoPhase), bic_2ph, a_2ph);
                    for i = 1:numel(acti_V)
                        acti_V{i}(twoPhase) = acti{i};
                    end
                end
            else
                [A_V, B_V, acti_V] = eos.getPhaseMixCoefficients(y, T, bic, a);
            end
            
            v{L_ix} = struct('acti', {acti_L}, 'A', {A_L}, 'B', {B_L}, 'a', {a}, 'bic', {bic}, 'Tr', {Tr});
            v{V_ix} = struct('acti', {acti_V}, 'A', {A_V}, 'B', {B_V}, 'a', {a}, 'bic', {bic}, 'Tr', {Tr});
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
