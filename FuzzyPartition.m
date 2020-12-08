classdef FuzzyPartition < handle
    
    properties(Constant)
        ALPHA = 0.5;
        ALPHA_AUX = sqrt(-2 * log(FuzzyPartition.ALPHA));
        I_SPD = 1;
        I_MU = 2;
    end
    
    properties
        age         (1,1) {mustBeReal, mustBeInteger, mustBeNonnegative} = 0;
        sim2Refit   (1,1) {mustBeReal, mustBePositive, mustBeLessThanOrEqual(sim2Refit,1)} = 0.7^(1/0.7^2);
        log2Smerge  (1,1) {mustBeReal, mustBeNonpositive} = log2(0.8)/0.8^2;
        maxMFs      (1,1) {mustBeReal, mustBeInteger, mustBePositive} = 25;
        mergedIDXs  (:,1) cell = cell(0,1); % {mustBeInteger, mustBePositive};
        mergedLgSim (:,1) {mustBeReal, mustBeNonpositive} = zeros(0,1);
        mergedMFs   (:,2) {mustBeReal} = zeros(0,2);
    end
    
    methods(Access = private)
        
        function log2Sim = log2PartitionSim(self)
            log2Sim = sum(self.mergedLgSim);
        end
        
        function sim = partitionSim(self)
            sim = 2^self.log2PartitionSim();
        end
        
        function len = totalMFs(self)
            len = length(self.mergedIDXs);
        end
        
        function bool = time2Refit(self, NC)
            bool = self.age > 5*(log2(NC)+1) ...
                || self.partitionSim() < self.sim2Refit;
        end
        
    end
    
    methods
        
        function self = FuzzyPartition(maxMFs, Smerge, sim2Refit, compMFs, weights)
            self.maxMFs = maxMFs;
            self.log2Smerge = log2(Smerge)/Smerge^2;
            self.sim2Refit = sim2Refit^(1/sim2Refit^2);
            components = self.calcAlphaSupport(compMFs, weights);
            self.refitPartition(components);
        end
        
        function updateSystem(self, compMFs, weights, thereIsANewComponent)
            components = FuzzyPartition.calcAlphaSupport(compMFs, weights);
            minAfter = Inf;
            for i_merged = length(self.mergedIDXs)-1:-1:1
                minAfter = min(minAfter, min(components.MF(self.mergedIDXs{i_merged+1}, self.I_MU)));
                if minAfter <= max(components.MF(self.mergedIDXs{i_merged}, self.I_MU))
                    self.refitPartition(components);
                    return;
                end
            end
            maxBefore = -Inf;
            for i_merged = 2:length(self.mergedIDXs)
                maxBefore = max(maxBefore, max(components.MF(self.mergedIDXs{i_merged-1}, self.I_MU)));
                if min(components.MF(self.mergedIDXs{i_merged}, self.I_MU)) <= maxBefore
                    self.refitPartition(components);
                    return;
                end
            end
            self.updateAllMergedMFs(components);
            self.updateAllmergedLgSim(components);
            self.age = self.age +1;
            if thereIsANewComponent
                self.age = self.age +1; % double-aging
                if self.time2Refit(length(components.weights))
                    self.refitPartition(components);
                    return;
                end
                self.fitNewComponent(components);
            end
            if self.time2Refit(length(components.weights))
                self.refitPartition(components);
            end
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
            for i_merged = 1:self.totalMFs()+1 % para cada MF da feature
                if i_merged <= self.totalMFs() && ...
                        self.mergedMFs(i_merged, self.I_MU) < newMF(self.I_MU)% se estiver antes do MF
                    continue;
                end
                if i_merged <= self.totalMFs() && ...
                        min(components.MF(self.mergedIDXs{i_merged}, self.I_MU)) <= newMF(self.I_MU)
                    idxNewMergedIDX = i_merged;
                    self.mergedIDXs{idxNewMergedIDX}(end+1,1) = idxNewMF;
                    self.updateIdxMergedMFs(components, idxNewMergedIDX);
                    self.updateIdxmergedLgSim(components, idxNewMergedIDX);
                elseif i_merged > 1 && ...
                        newMF(self.I_MU) <= max(components.MF(self.mergedIDXs{i_merged-1}, self.I_MU))
                    idxNewMergedIDX = i_merged-1;
                    self.mergedIDXs{idxNewMergedIDX}(end+1,1) = idxNewMF;
                    self.updateIdxMergedMFs(components, idxNewMergedIDX);
                    self.updateIdxmergedLgSim(components, idxNewMergedIDX);
                else % it is between
                    self.putNewMergedIDXAt({idxNewMF}, i_merged, components);
                    self.tryMergeMFs(components);
                end
                break;
            end
        end
        
        function tryMergeMFs(self, components)
            if self.totalMFs() == 1
                return;
            end
            log2SimOnMerge = zeros(length(self.mergedIDXs)-1, 1);
            for ii = 1:length(log2SimOnMerge) % para cada par de componentes vizinhos
                log2SimOnMerge(ii) = self.log2SimilarityOnMergeIdx([ii, ii+1], components);
            end
            log2SimImprov = log2SimOnMerge - (self.mergedLgSim(1:end-1) + self.mergedLgSim(2:end));
            [improv, idxMerged] = max(log2SimImprov);
            while ~isempty(log2SimImprov) && ...
                    (   improv >= 0 || self.totalMFs() > self.maxMFs ...
                    ||  self.log2PartitionSim() + improv > self.log2Smerge )
                self.mergedIDXs{idxMerged} = [self.mergedIDXs{idxMerged}; self.mergedIDXs{idxMerged+1}];
                self.popMergedPropsAt(idxMerged+1);
                self.updateIdxMergedMFs(components, idxMerged);
                self.mergedLgSim(idxMerged) = log2SimOnMerge(idxMerged);
                log2SimOnMerge(idxMerged) = [];
                log2SimImprov(idxMerged) = [];
                if idxMerged > 1
                    log2SimOnMerge(idxMerged-1) = self.log2SimilarityOnMergeIdx([idxMerged-1, idxMerged], components);
                    log2SimImprov(idxMerged-1) = log2SimOnMerge(idxMerged-1) ...
                        - (self.mergedLgSim(idxMerged-1) + self.mergedLgSim(idxMerged));
                end
                if idxMerged < length(self.mergedIDXs)
                    log2SimOnMerge(idxMerged) = self.log2SimilarityOnMergeIdx([idxMerged, idxMerged+1], components);
                    log2SimImprov(idxMerged) = log2SimOnMerge(idxMerged) ...
                        - (self.mergedLgSim(idxMerged) + self.mergedLgSim(idxMerged+1));
                end
                [improv, idxMerged] = max(log2SimImprov);
            end
        end
        
        function refitPartition(self, components)
            [~, sortedMu] = sort(components.MF(:, self.I_MU));
            self.mergedIDXs = num2cell(sortedMu);
            self.updateAllMergedMFs(components);
            self.updateAllmergedLgSim(components);
            self.tryMergeMFs(components);
            self.age = 0;
        end
        
        function putNewMergedIDXAt(self, newMergedIDX, idxNewMergedIDX, components)
            self.mergedIDXs = [ self.mergedIDXs(1:idxNewMergedIDX-1); ...
                newMergedIDX;   self.mergedIDXs(idxNewMergedIDX:end) ];
            self.mergedMFs(idxNewMergedIDX+1:end+1, :) = self.mergedMFs(idxNewMergedIDX:end, :);
            self.mergedLgSim(idxNewMergedIDX+1:end+1, :) = self.mergedLgSim(idxNewMergedIDX:end, :);
            adjIdxNewMergedIDX = max(idxNewMergedIDX-1, 1):min(self.totalMFs(), idxNewMergedIDX+1);
            self.updateIdxMergedMFs(components, adjIdxNewMergedIDX);
            self.updateIdxmergedLgSim(components, adjIdxNewMergedIDX);
        end
        
        function popMergedPropsAt(self, idxMerged)
            self.mergedIDXs(idxMerged) = [];
            self.mergedLgSim(idxMerged) = [];
            self.mergedMFs(idxMerged,:) = [];
        end
        
        function updateAllmergedLgSim(self, components)
            self.mergedLgSim = zeros(size(self.mergedMFs,1), 1);
            self.updateIdxmergedLgSim(components, 1:size(self.mergedMFs,1));
        end
        
        function updateIdxmergedLgSim(self, components, idxMergeds)
            for idx = idxMergeds
                sim = self.similarity( self.mergedMFs(idx,:), ...
                    components.MF(self.mergedIDXs{idx}, :) );
                weight = components.weights(self.mergedIDXs{idx});
                self.mergedLgSim(idx) = sum(weight .* log2(sim) ./ sim.^2);
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
            merged = FuzzyPartition.mergeAlphaSupport( [ minXroot maxXroot ] );
            sim = self.similarity( merged, components.MF(idxs, :) );
            weight = components.weights(idxs);
            log2sim = sum(weight .* log2(sim) ./ sim.^2);
        end
        
        function sim = similarity(self, mergedOne, subMFs)
            mudiff = abs(mergedOne(self.I_MU) - subMFs(:, self.I_MU));
            spdaux = [ mergedOne(self.I_SPD) + zeros(size(subMFs, 1), 1), ...
                subMFs(:, self.I_SPD) ];
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
        
    end
    
    methods(Static, Access = private)
        
        function components = calcAlphaSupport(compMFs, weights)
            components.MF = compMFs;
            components.weights = weights(:);
            aux = components.MF(:, FuzzyPartition.I_SPD) * FuzzyPartition.ALPHA_AUX;
            components.alphaSupport = [ ... 
                ( components.MF(:, FuzzyPartition.I_MU) - aux ), ...
                ( components.MF(:, FuzzyPartition.I_MU) + aux ) ];
        end
        
        function merged = mergeAlphaSupport(alphaSupports)
            minXroot = min(alphaSupports(:, 1));
            maxXroot = max(alphaSupports(:, 2));
            newMu = (minXroot + maxXroot) / 2;
            newSigma = (maxXroot - newMu) / FuzzyPartition.ALPHA_AUX;
            newSigma = max(newSigma, 1e-7);
            merged = [newSigma, newMu];
        end
        
    end
end

