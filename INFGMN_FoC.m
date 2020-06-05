classdef INFGMN_FoC < handle
    
    properties(Constant)
        ALPHA = 0.5;
        ALPHA_AUX = sqrt(-2 * log(INFGMN_FoC.ALPHA));
        I_SPD = 1;
        I_MU = 2;
    end
    
    properties
        Sage    (1,1) {mustBeReal, mustBeNonnegative, mustBeLessThanOrEqual(Sage,1)} = 1;
        Smerge  (1,1) {mustBeReal, mustBePositive, mustBeLessThan(Smerge,1)} = 0.4;
        vmax    (1,1) {mustBeReal, mustBePositive, mustBeLessThan(vmax,1)} = 0.8;
        maxFoCSize  (1,1) {mustBeReal, mustBeInteger, mustBePositive} = 7;
%         vweightexp  (1,1) {mustBeReal} = 0;
        mergedIDXs  (:,1) cell = cell(0,1); %{mustBeInteger, mustBePositive};
        mergedSim   (:,1) {mustBeReal, mustBeNonnegative, mustBeLessThanOrEqual(mergedSim,1)} = zeros(0,1);
        mergedMFs   (:,2) {mustBeReal} = zeros(0,2);
    end
    
    methods
        
        function self = INFGMN_FoC(maxFoCSize, Smerge, vmax, compMFs)
            self.maxFoCSize = maxFoCSize;
            self.Smerge = Smerge;
            self.vmax = vmax;
%             vmerge = (1+vmax)/2;
%             self.vweightexp = log(Smerge)/log(vmerge)-1;
            components = self.calcAlphaSupport(compMFs);
            self.renewFoC(components);
        end
        
        function sim = FoCSim(self)
            sim = prod(self.mergedSim) ^ (1/ length(vertcat(self.mergedIDXs{:})) );
        end
        
        function len = FoCSize(self)
            len = length(self.mergedIDXs);
        end
        
    end
    
    methods
        
        function updateSystem(self, compMFs, thereIsANewComponent)
            components = INFGMN_FoC.calcAlphaSupport(compMFs);
            minAfter = Inf;
            for i_merged = length(self.mergedIDXs)-1:-1:1
                minAfter = min(minAfter, min(components.MF(self.mergedIDXs{i_merged+1}, self.I_MU)));
                if minAfter <= max(components.MF(self.mergedIDXs{i_merged}, self.I_MU))
                    self.renewFoC(components);
                    return;
                end
            end
            maxBefore = -Inf;
            for i_merged = 2:length(self.mergedIDXs)
                maxBefore = max(maxBefore, max(components.MF(self.mergedIDXs{i_merged-1}, self.I_MU)));
                if min(components.MF(self.mergedIDXs{i_merged}, self.I_MU)) <= maxBefore
                    self.renewFoC(components);
                    return;
                end
            end
            self.updateAllMergedMFs(components);
            self.updateAllMergedSim(components);
            self.Sage = self.aging();
            if thereIsANewComponent
                if self.aging() < self.Smerge
                    self.renewFoC(components);
                    return;
                end
                self.fitNewComponent(components);
                self.Sage = self.aging(); % double-aging
            end
            if self.Sage < self.Smerge
                self.renewFoC(components);
            end
        end
        
        function age = aging(self)
            age = sqrt(min(self.Sage, self.Smerge)) * self.FoCSim();
        end
        
        function purge(self, idx)
            for idxMerged = length(self.mergedIDXs):-1:1
                self.mergedIDXs{idxMerged} = setdiff(self.mergedIDXs{idxMerged}, idx);
                if isempty(self.mergedIDXs{idxMerged})
                    self.popMergedPropsAt(idxMerged);
                else
                    filtro = self.mergedIDXs{idxMerged} > idx;
                    self.mergedIDXs{idxMerged}(filtro) = ...
                        self.mergedIDXs{idxMerged}(filtro) - 1;
                end
            end
        end
        
    end
    
    methods(Access = private)
        
        function fitNewComponent(self, components)
            newMF = components.MF(end, :);
            idxNewMF = size(components.MF, 1);
            for i_merged = 1:self.FoCSize()+1 % para cada MF da feature
                if i_merged <= self.FoCSize() && ...
                        self.mergedMFs(i_merged, self.I_MU) < newMF(self.I_MU)% se estiver antes do MF
                    continue;
                end
                if i_merged <= self.FoCSize() && ...
                        min(components.MF(self.mergedIDXs{i_merged}, self.I_MU)) <= newMF(self.I_MU)
                    idxNewMergedIDX = i_merged;
                    self.mergedIDXs{idxNewMergedIDX}(end+1,1) = idxNewMF;
                    self.updateIdxMergedMFs(components, idxNewMergedIDX);
                    self.updateIdxMergedSim(components, idxNewMergedIDX);
                elseif i_merged > 1 && ...
                        newMF(self.I_MU) <= max(components.MF(self.mergedIDXs{i_merged-1}, self.I_MU))
                    idxNewMergedIDX = i_merged-1;
                    self.mergedIDXs{idxNewMergedIDX}(end+1,1) = idxNewMF;
                    self.updateIdxMergedMFs(components, idxNewMergedIDX);
                    self.updateIdxMergedSim(components, idxNewMergedIDX);
                else % it is between
                    self.putNewMergedIDXAt({idxNewMF}, i_merged, components);
                    self.tryMergeMFs(components);
                end
                break;
            end
        end
        
        function tryMergeMFs(self, components)
            if self.FoCSize() == 1
                return;
            end
            simOnMerge = zeros(length(self.mergedIDXs)-1, 1);
            for ii = 1:length(simOnMerge) % para cada par de componentes vizinhos
                simOnMerge(ii) = self.similarIdx([ii, ii+1], components);
            end
            simImprov = simOnMerge ./ (self.mergedSim(1:end-1) .* self.mergedSim(2:end));
            [improv, idxMerged] = max(simImprov);
            while ~isempty(simImprov) && ( improv >= 1 || self.FoCSize() > self.maxFoCSize )
                self.mergedIDXs{idxMerged} = [self.mergedIDXs{idxMerged}; self.mergedIDXs{idxMerged+1}];
                self.popMergedPropsAt(idxMerged+1);
                self.updateIdxMergedMFs(components, idxMerged);
                self.mergedSim(idxMerged) = simOnMerge(idxMerged);
                simOnMerge(idxMerged) = [];
                simImprov(idxMerged) = [];
                if idxMerged > 1
                    simOnMerge(idxMerged-1) = self.similarIdx([idxMerged-1, idxMerged], components);
                    simImprov(idxMerged-1) = simOnMerge(idxMerged-1) ...
                        / (self.mergedSim(idxMerged-1) * self.mergedSim(idxMerged));
                end
                if idxMerged < length(self.mergedIDXs)
                    simOnMerge(idxMerged) = self.similarIdx([idxMerged, idxMerged+1], components);
                    simImprov(idxMerged) = simOnMerge(idxMerged) ...
                        / (self.mergedSim(idxMerged) * self.mergedSim(idxMerged+1));
                end
                [improv, idxMerged] = max(simImprov);
            end
        end
        
        function renewFoC(self, components)
            [~, sortedMu] = sort(components.MF(:, self.I_MU));
            self.mergedIDXs = num2cell(sortedMu);
            self.updateAllMergedMFs(components);
            self.updateAllMergedSim(components);
            self.tryMergeMFs(components);
            self.Sage = 1;
        end
        
        function putNewMergedIDXAt(self, newMergedIDX, idxNewMergedIDX, components)
            self.mergedIDXs = [ self.mergedIDXs(1:idxNewMergedIDX-1); ...
                newMergedIDX;   self.mergedIDXs(idxNewMergedIDX:end) ];
            self.mergedMFs(idxNewMergedIDX+1:end+1, :) = self.mergedMFs(idxNewMergedIDX:end, :);
            self.mergedSim(idxNewMergedIDX+1:end+1, :) = self.mergedSim(idxNewMergedIDX:end, :);
            adjIdxNewMergedIDX = max(idxNewMergedIDX-1, 1):min(self.FoCSize(), idxNewMergedIDX+1);
            self.updateIdxMergedMFs(components, adjIdxNewMergedIDX);
            self.updateIdxMergedSim(components, adjIdxNewMergedIDX);
        end
        
        function popMergedPropsAt(self, idxMerged)
            self.mergedIDXs(idxMerged) = [];
            self.mergedSim(idxMerged) = [];
            self.mergedMFs(idxMerged,:) = [];
        end
        
        function updateAllMergedSim(self, components)
            self.mergedSim = zeros(size(self.mergedMFs,1), 1);
            self.updateIdxMergedSim(components, 1:size(self.mergedMFs,1));
        end
        
        function updateIdxMergedSim(self, components, idxMergeds)
            for idx = idxMergeds
                sim = self.similarity( self.mergedMFs(idx,:), ...
                    components.MF(self.mergedIDXs{idx}, :) );
                self.mergedSim(idx) = prod(sim);
            end
        end
        
        function updateAllMergedMFs(self, components)
            self.mergedMFs = zeros(length(self.mergedIDXs), 2);
            self.updateIdxMergedMFs(components, 1:length(self.mergedIDXs));
        end
        
        function updateIdxMergedMFs(self, components, idxMergeds)
            for idx = idxMergeds
                minXroot = min(components.alphaSupport( self.mergedIDXs{idx}, 1 ));
                if idx > 1
                    minXroot = max( minXroot, ...
                        max(components.MF(self.mergedIDXs{idx-1}, self.I_MU)) );
                end
                maxXroot = max(components.alphaSupport( self.mergedIDXs{idx}, 2 ));
                if idx < length(self.mergedIDXs)
                    maxXroot = min( maxXroot, ...
                        min(components.MF(self.mergedIDXs{idx+1}, self.I_MU)) );
                end
                self.mergedMFs(idx, :) = self.mergeAlphaSupport( [ minXroot maxXroot ] );
            end
        end
        
        function [sim, merged] = similarIdx(self, idxMergeds, components)
            idxs = vertcat(self.mergedIDXs{idxMergeds});
            minXroot = min(components.alphaSupport( idxs, 1 ));
            if min(idxMergeds) > 1
                minXroot = max( minXroot, ...
                    max(components.MF(self.mergedIDXs{min(idxMergeds)-1}, self.I_MU)) );
            end
            maxXroot = max(components.alphaSupport( idxs, 2 ));
            if max(idxMergeds) < length(self.mergedIDXs)
                maxXroot = min( maxXroot, ...
                    min(components.MF(self.mergedIDXs{max(idxMergeds)+1}, self.I_MU)) );
            end
            merged = INFGMN_FoC.mergeAlphaSupport( [ minXroot maxXroot ] );
            sim = self.similarity( merged, components.MF(idxs, :) );
            sim = prod(sim);
        end
        
        function sim = similarity(self, mergedOne, subMFs)
            mudiff = abs(mergedOne(self.I_MU) - subMFs(:, self.I_MU));
            spdaux = [ mergedOne(self.I_SPD) + zeros(size(subMFs, 1), 1), ...
                subMFs(:, self.I_SPD) ];
%             spdaux = [ zeros(size(subMFs, 1), 1) subMFs(:, self.I_SPD) ];
%             for ii = 1:size(subMFs,1)
%                 spdaux(ii, 1) = min( ...
%                     mergedOne(self.I_SPD), ...
%                     A(self.I_MU) - A(self.I_SPD) * self.ALPHA_AUX ...
%                     );
%             end
            spdmin = min( spdaux, [], 2 );
            spdmax = max( spdaux, [], 2 );
            spddiff = spdmax - spdmin;
            spdsum  = spdmax + spdmin;
            sqrtneglnw = mudiff ./ spddiff;
            sqrtneglnv = mudiff ./ spdsum; % v = possibility
            ohm = spddiff .* erf(sqrtneglnw) ...
                - spdsum  .* erf(sqrtneglnv);
            sim  = (2 * spdmin + ohm) ...
                ./ (2 * spdmax - ohm);
            sim(mudiff==0 & spddiff==0) = 1;
            assert(~any(isnan(sim)));
            lim_inf = 1e-15;
            sim(-1e-15 < sim & sim < lim_inf) = lim_inf;
            assert(all(sim > 0));
        end
        
%         function [sim, v] = similarity(self, A, B)
%             v = self.possibility(A, B);
%             vWeight = v^self.vweightexp;
%             if vWeight * v >= self.Smerge
%                 sim = v;
%             else
%                 spdmin = min(A(self.I_SPD), B(self.I_SPD));
%                 beta = 2 * spdmin / (A(self.I_SPD) + B(self.I_SPD));
%                 sqrtneglnv = sqrt(-log(v));
%                 psi = beta ...
%                     + (1 - beta) * erf( (1/(1-beta)) * sqrtneglnv ) ...
%                     - erf( sqrtneglnv );
%                 sim = vWeight * v ...
%                     + (1 - vWeight) * (psi / (2 - psi));
%             end
%         end
        
    end
    
    methods(Static, Access = private)
        
        function components = calcAlphaSupport(compMFs)
            components.MF = compMFs;
            aux = components.MF(:, INFGMN_FoC.I_SPD) * INFGMN_FoC.ALPHA_AUX;
            components.alphaSupport = [ ... 
                ( components.MF(:, INFGMN_FoC.I_MU) - aux ), ...
                ( components.MF(:, INFGMN_FoC.I_MU) + aux ) ];
        end
        
        function merged = mergeAlphaSupport(alphaSupports)
            minXroot = min(alphaSupports(:, 1));
            maxXroot = max(alphaSupports(:, 2));
            newMu = (minXroot + maxXroot) / 2;
            newSigma = (maxXroot - newMu) / INFGMN_FoC.ALPHA_AUX;
            newSigma = max(newSigma, 1e-7);
            merged = [newSigma, newMu];
        end
        
%         function v = possibility(A, B)
%             v = exp( ...
%                 -( (A(INFGMN_FoC.I_MU)  - B(INFGMN_FoC.I_MU)) ...
%                 /  (A(INFGMN_FoC.I_SPD) + B(INFGMN_FoC.I_SPD)) ...
%                 )^2 /2 );
%         end
        
    end
end

