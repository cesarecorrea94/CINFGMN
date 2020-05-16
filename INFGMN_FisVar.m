classdef INFGMN_FisVar < handle
    
    properties(Constant)
        ALPHA = 0.5;
        ALPHA_AUX = sqrt(-2 * log(INFGMN_FisVar.ALPHA));
        I_SPD = 1;
        I_MU = 2;
    end
    
    properties
        Smerge  (1,1) {mustBeReal, mustBePositive, mustBeLessThan(Smerge,1)} = 0.4;
        vmax    (1,1) {mustBeReal, mustBePositive, mustBeLessThan(vmax,1)} = 0.8;
        vweightexp  (1,1) {mustBeReal} = 0;
        mergedIDXs  (:,1) cell = cell(0,1); %{mustBeInteger, mustBePositive};
        mergedMFs   (:,2) {mustBeReal} = zeros(0,2);
    end
    
    methods
        
        function self = INFGMN_FisVar(Smerge, vmax, compMFs)
            self.Smerge = Smerge;
            self.vmax = vmax;
            vmerge = (1+vmax)/2;
            self.vweightexp = log(Smerge)/log(vmerge)-1;
            [~, sortedMu] = sort(compMFs(:, self.I_MU));
            self.mergedIDXs = num2cell(sortedMu);
            self.mergedMFs = compMFs(sortedMu, :);
            self.tryMergeMFs(compMFs, self.calcAlphaSupport(compMFs));
        end
        
        function fitNewComponent(self, compMFs)
            alphaSupport = self.updateAllMergedMFs(compMFs);
            newMF = compMFs(end, :);
            idxNewMF = size(compMFs, 1);
            for i_merged = 1:length(self.mergedIDXs)+1 % para cada MF da feature
                if i_merged <= length(self.mergedIDXs) && ...
                        self.mergedMFs(i_merged, self.I_MU) < newMF(self.I_MU)% se estiver antes do MF
                    continue;
                end
                if i_merged > length(self.mergedIDXs)
                    idxNewMergedIDX = length(self.mergedIDXs);
                elseif i_merged == 1
                    idxNewMergedIDX = 1;
                elseif newMF(self.I_MU) <= max(compMFs(self.mergedIDXs{i_merged-1}, self.I_MU))
                    idxNewMergedIDX = i_merged-1;
                elseif min(compMFs(self.mergedIDXs{i_merged}, self.I_MU)) <= newMF(self.I_MU)
                    idxNewMergedIDX = i_merged;
                else % it is between
                    pre_sim = self.similarIdx(i_merged-1, idxNewMF, compMFs);
                    pos_sim = self.similarIdx(i_merged, idxNewMF, compMFs);
                    if any([pre_sim, pos_sim] > self.Smerge)
                        if pre_sim > pos_sim
                            idxNewMergedIDX = i_merged-1;
                        else
                            idxNewMergedIDX = i_merged;
                        end
                    else % merge both to disintegrate soon after
                        self.mergedIDXs{i_merged-1} = [ ...
                            self.mergedIDXs{i_merged-1}; ...
                            self.mergedIDXs{i_merged} ];
                        self.mergedIDXs(i_merged) = [];
                        self.mergedMFs(i_merged, :) = [];
                        idxNewMergedIDX = i_merged-1;
                    end
                end
                self.mergedIDXs{idxNewMergedIDX}(end+1,1) = idxNewMF;
                self.updateIdxMergedMFs(alphaSupport, idxNewMergedIDX);
                % check if is similar.
                break;
            end
        end
        
        function updateSystem(self, compMFs)
            self.updateAllMergedMFs(compMFs);
            disintegrate = false(length(self.mergedIDXs), 1);
            minAfter = Inf;
            for i_merged = length(self.mergedIDXs)-1:-1:1
                minAfter = min(minAfter, min(compMFs(self.mergedIDXs{i_merged+1}, self.I_MU)));
                if minAfter <= max(compMFs(self.mergedIDXs{i_merged}, self.I_MU))
                    disintegrate(i_merged) = true;
                end
            end
            maxBefore = -Inf;
            for i_merged = 2:length(self.mergedIDXs)
                maxBefore = max(maxBefore, max(compMFs(self.mergedIDXs{i_merged-1}, self.I_MU)));
                if min(compMFs(self.mergedIDXs{i_merged}, self.I_MU)) <= maxBefore
                    disintegrate(i_merged) = true;
                end
            end
            for i_merged = 1:length(self.mergedIDXs)
                if ~disintegrate(i_merged) && ...
                        self.similarIdx( i_merged, [], compMFs ) < self.Smerge
                    disintegrate(i_merged) = true;
                end
            end
            if any(disintegrate)
                self.mergedIDXs = [self.mergedIDXs(~disintegrate); ...
                    num2cell(vertcat(self.mergedIDXs{disintegrate})) ];
                alphaSupport = self.updateAllMergedMFs(compMFs);
                [~, sortedMu] = sort(self.mergedMFs(:, self.I_MU));
                self.mergedIDXs = self.mergedIDXs(sortedMu);
                self.mergedMFs = self.mergedMFs(sortedMu, :);
                self.tryMergeMFs(compMFs, alphaSupport);
            end
        end
        
        function purge(self, idx)
            for idxMerged = length(self.mergedIDXs):-1:1
                self.mergedIDXs{idxMerged} = setdiff(self.mergedIDXs{idxMerged}, idx);
                if isempty(self.mergedIDXs{idxMerged})
                    self.mergedIDXs(idxMerged) = [];
                else
                    filtro = self.mergedIDXs{idxMerged} > idx;
                    self.mergedIDXs{idxMerged}(filtro) = ...
                        self.mergedIDXs{idxMerged}(filtro) - 1;
                end
            end
        end
        
    end
    
    methods(Access = private)
        
        function alphaSupport = updateAllMergedMFs(self, compMFs)
            alphaSupport = self.calcAlphaSupport(compMFs);
            self.mergedMFs = zeros(length(self.mergedIDXs), 2);
            self.updateIdxMergedMFs(alphaSupport, 1:length(self.mergedIDXs));
        end
        
        function updateIdxMergedMFs(self, alphaSupport, idxs)
            for ii = idxs
                self.mergedMFs(ii, :) = self.mergeAlphaSupport( alphaSupport( self.mergedIDXs{ii}, : ) );
            end
        end
        
        function tryMergeMFs(self, compMFs, alphaSupport)
            similarities = Inf(length(self.mergedIDXs)-1, 1);
            for ii = 1:length(similarities) % para cada par de componentes vizinhos
                similarities(ii) = self.similarIdx([ii, ii+1], [], compMFs);
            end
            [sim, idxMerged] = max(similarities);
            while sim >= self.Smerge
                self.mergedIDXs{idxMerged} = [self.mergedIDXs{idxMerged}; self.mergedIDXs{idxMerged+1}];
                self.updateIdxMergedMFs(alphaSupport, idxMerged);
                self.mergedIDXs(idxMerged+1) = [];
                self.mergedMFs(idxMerged+1, :) = [];
                similarities(idxMerged) = [];
                if idxMerged > 1
                    similarities(idxMerged-1) = self.similarIdx([idxMerged-1, idxMerged], [], compMFs);
                end
                if idxMerged < length(self.mergedIDXs)
                    similarities(idxMerged) = self.similarIdx([idxMerged, idxMerged+1], [], compMFs);
                end
                [sim, idxMerged] = max(similarities);
            end
        end
        
        function [sim, v] = similarIdx(self, idxMergeds, idxMFs, compMFs)
            idxs = [ vertcat(self.mergedIDXs{idxMergeds}); idxMFs(:) ];
            merged = INFGMN_FisVar.mergeAlphaSupport( ...
                INFGMN_FisVar.calcAlphaSupport( ...
                [ self.mergedMFs(idxMergeds, :); compMFs(idxMFs, :) ] ));
            sqrtneglnv = abs( merged(INFGMN_FisVar.I_MU)  - compMFs(idxs, INFGMN_FisVar.I_MU) ) ...
                ./          ( merged(INFGMN_FisVar.I_SPD) + compMFs(idxs, INFGMN_FisVar.I_SPD) );
            [~, idxMinV] = max(sqrtneglnv);
            [sim, v] = self.similarity( merged, compMFs(idxs(idxMinV), :) );
        end
        
        function [sim, v] = similarity(self, A, B)
            v = self.possibility(A, B);
            vWeight = v^self.vweightexp;
            if vWeight * v >= self.Smerge
                sim = v;
            else
                spdmin = min(A(INFGMN_FisVar.I_SPD), B(INFGMN_FisVar.I_SPD));
                beta = 2 * spdmin / (A(INFGMN_FisVar.I_SPD) + B(INFGMN_FisVar.I_SPD));
                sqrtneglnv = sqrt(-log(v));
                psi = beta ...
                    + (1 - beta) * erf( (1/(1-beta)) * sqrtneglnv ) ...
                    - erf( sqrtneglnv );
                sim = vWeight * v ...
                    + (1 - vWeight) * (psi / (2 - psi));
            end
        end
        
    end
    
    methods(Static, Access = private)
        
        function alphaSupport = calcAlphaSupport(compMFs)
            aux = compMFs(:, INFGMN_FisVar.I_SPD) * INFGMN_FisVar.ALPHA_AUX;
            alphaSupport = [ ...
                ( compMFs(:, INFGMN_FisVar.I_MU) - aux ), ...
                ( compMFs(:, INFGMN_FisVar.I_MU) + aux ) ];
        end
        
        function merged = mergeAlphaSupport(alphaSupports)
            minXroot = min(alphaSupports(:, 1));
            maxXroot = max(alphaSupports(:, 2));
            newMu = (minXroot + maxXroot) / 2;
            newSigma = (maxXroot - newMu) / INFGMN_FisVar.ALPHA_AUX;
            merged = [newSigma, newMu];
        end
        
        function v = possibility(A, B)
            v = exp( ...
                -( (A(INFGMN_FisVar.I_MU)  - B(INFGMN_FisVar.I_MU)) ...
                /  (A(INFGMN_FisVar.I_SPD) + B(INFGMN_FisVar.I_SPD)) ...
                )^2 /2 );
        end
        
    end
end

