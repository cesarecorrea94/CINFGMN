classdef INFGMN_FoC < handle
    
    properties(Constant)
        ALPHA = 0.5;
        ALPHA_AUX = sqrt(-2 * log(INFGMN_FoC.ALPHA));
        I_SPD = 1;
        I_MU = 2;
    end
    
    properties
        log2Sage    (1,1) {mustBeReal, mustBeNonpositive} = 0;
        log2Smerge  (1,1) {mustBeReal, mustBeNonpositive} = log2(0.85);
        log2Sdeath  (1,1) {mustBeReal, mustBeNonpositive} = log2(0.7);
        vmax    (1,1) {mustBeReal, mustBePositive, mustBeLessThan(vmax,1)} = 0.8;
        maxFoCSize  (1,1) {mustBeReal, mustBeInteger, mustBePositive} = 7;
%         vweightexp  (1,1) {mustBeReal} = 0;
        mergedIDXs      (:,1) cell = cell(0,1); %{mustBeInteger, mustBePositive};
        mergedLog2Sim   (:,1) {mustBeReal, mustBeNonpositive} = zeros(0,1);
        mergedMFs       (:,2) {mustBeReal} = zeros(0,2);
    end
    
    methods
        
        function self = INFGMN_FoC(maxFoCSize, Smerge, Sdeath, vmax, compMFs, weights)
            self.maxFoCSize = maxFoCSize;
            self.log2Smerge = log2(Smerge)/Smerge^2;
            self.log2Sdeath = log2(Sdeath)/Sdeath^2;
            self.vmax = vmax;
%             vmerge = (1+vmax)/2;
%             self.vweightexp = log(Smerge)/log(vmerge)-1;
            components = self.calcAlphaSupport(compMFs, weights);
            self.renewFoC(components);
        end
        
        function log2Sim = log2FoCSim(self)
            log2Sim = sum(self.mergedLog2Sim);
        end
        
        function len = FoCSize(self)
            len = length(self.mergedIDXs);
        end
        
    end
    
    methods
        
        function updateSystem(self, compMFs, weights, thereIsANewComponent)
            components = INFGMN_FoC.calcAlphaSupport(compMFs, weights);
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
            self.updateAllMergedLog2Sim(components);
            self.log2Sage = self.log2Aging();
            if thereIsANewComponent
                if self.log2Aging() < self.log2Sdeath
                    self.renewFoC(components);
                    return;
                end
                self.fitNewComponent(components);
                self.log2Sage = self.log2Aging(); % double-aging
            end
            if self.log2Sage < self.log2Sdeath
                self.renewFoC(components);
            end
        end
        
        function log2SAge = log2Aging(self)
            log2SAge = self.log2Sage/2 + min(self.log2FoCSim(), self.log2Sdeath/2);
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
                    self.updateIdxMergedLog2Sim(components, idxNewMergedIDX);
                elseif i_merged > 1 && ...
                        newMF(self.I_MU) <= max(components.MF(self.mergedIDXs{i_merged-1}, self.I_MU))
                    idxNewMergedIDX = i_merged-1;
                    self.mergedIDXs{idxNewMergedIDX}(end+1,1) = idxNewMF;
                    self.updateIdxMergedMFs(components, idxNewMergedIDX);
                    self.updateIdxMergedLog2Sim(components, idxNewMergedIDX);
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
            log2SimOnMerge = zeros(length(self.mergedIDXs)-1, 1);
            for ii = 1:length(log2SimOnMerge) % para cada par de componentes vizinhos
                log2SimOnMerge(ii) = self.log2SimilarityOnMergeIdx([ii, ii+1], components);
            end
            log2SimImprov = log2SimOnMerge - (self.mergedLog2Sim(1:end-1) + self.mergedLog2Sim(2:end));
            [improv, idxMerged] = max(log2SimImprov);
            while ~isempty(log2SimImprov) && ...
                    (   improv >= 0 || self.FoCSize() > self.maxFoCSize ...
                    ||  self.log2FoCSim() + improv > self.log2Smerge )
                self.mergedIDXs{idxMerged} = [self.mergedIDXs{idxMerged}; self.mergedIDXs{idxMerged+1}];
                self.popMergedPropsAt(idxMerged+1);
                self.updateIdxMergedMFs(components, idxMerged);
                self.mergedLog2Sim(idxMerged) = log2SimOnMerge(idxMerged);
                log2SimOnMerge(idxMerged) = [];
                log2SimImprov(idxMerged) = [];
                if idxMerged > 1
                    log2SimOnMerge(idxMerged-1) = self.log2SimilarityOnMergeIdx([idxMerged-1, idxMerged], components);
                    log2SimImprov(idxMerged-1) = log2SimOnMerge(idxMerged-1) ...
                        - (self.mergedLog2Sim(idxMerged-1) + self.mergedLog2Sim(idxMerged));
                end
                if idxMerged < length(self.mergedIDXs)
                    log2SimOnMerge(idxMerged) = self.log2SimilarityOnMergeIdx([idxMerged, idxMerged+1], components);
                    log2SimImprov(idxMerged) = log2SimOnMerge(idxMerged) ...
                        - (self.mergedLog2Sim(idxMerged) + self.mergedLog2Sim(idxMerged+1));
                end
                [improv, idxMerged] = max(log2SimImprov);
            end
        end
        
        function renewFoC(self, components)
            [~, sortedMu] = sort(components.MF(:, self.I_MU));
            self.mergedIDXs = num2cell(sortedMu);
            self.updateAllMergedMFs(components);
            self.updateAllMergedLog2Sim(components);
            self.tryMergeMFs(components);
            self.log2Sage = 0;
        end
        
        function putNewMergedIDXAt(self, newMergedIDX, idxNewMergedIDX, components)
            self.mergedIDXs = [ self.mergedIDXs(1:idxNewMergedIDX-1); ...
                newMergedIDX;   self.mergedIDXs(idxNewMergedIDX:end) ];
            self.mergedMFs(idxNewMergedIDX+1:end+1, :) = self.mergedMFs(idxNewMergedIDX:end, :);
            self.mergedLog2Sim(idxNewMergedIDX+1:end+1, :) = self.mergedLog2Sim(idxNewMergedIDX:end, :);
            adjIdxNewMergedIDX = max(idxNewMergedIDX-1, 1):min(self.FoCSize(), idxNewMergedIDX+1);
            self.updateIdxMergedMFs(components, adjIdxNewMergedIDX);
            self.updateIdxMergedLog2Sim(components, adjIdxNewMergedIDX);
        end
        
        function popMergedPropsAt(self, idxMerged)
            self.mergedIDXs(idxMerged) = [];
            self.mergedLog2Sim(idxMerged) = [];
            self.mergedMFs(idxMerged,:) = [];
        end
        
        function updateAllMergedLog2Sim(self, components)
            self.mergedLog2Sim = zeros(size(self.mergedMFs,1), 1);
            self.updateIdxMergedLog2Sim(components, 1:size(self.mergedMFs,1));
        end
        
        function updateIdxMergedLog2Sim(self, components, idxMergeds)
            for idx = idxMergeds
                sim = self.similarity( self.mergedMFs(idx,:), ...
                    components.MF(self.mergedIDXs{idx}, :) );
                weight = components.weights(self.mergedIDXs{idx});
                self.mergedLog2Sim(idx) = sum(weight .* log2(sim) ./ sim.^2);
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
        
        function [log2sim, merged] = log2SimilarityOnMergeIdx(self, idxMergeds, components)
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
            weight = components.weights(idxs);
            log2sim = sum(weight .* log2(sim) ./ sim.^2);
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
        
        function components = calcAlphaSupport(compMFs, weights)
            components.MF = compMFs;
            components.weights = weights(:);
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

