classdef INFGMN_FisVar < handle
    
    properties(Constant)
        ALPHA = 0.5;
        ALPHA_AUX = sqrt(-2 * log(INFGMN_FisVar.ALPHA));
        I_SPD = 1;
        I_MU = 2;
    end
    
    properties
        Smerge;
        vmax;
        mergedIDXs;
        mergedMFs;
    end
    
    methods
        
        function self = INFGMN_FisVar(Smerge, vmax, compsMFs)
            self.Smerge = Smerge;
            self.vmax = vmax;
            self.renewMFs(compsMFs);
        end
        
        function fitNewComponent(self, compMFs)
            alphaSupport = self.updateMergedMFs(compMFs);
            newMF = compMFs(end, :);
            idxNewMF = size(compMFs, 1);
            if  newMF(self.I_MU) <= self.mergedMFs(1, self.I_MU) % se estiver antes do primeiro MF
                if self.similarIdx( 1, idxNewMF, compMFs ) >= self.Smerge
                    self.mergedIDXs{1}(end+1) = idxNewMF;
                    self.updateMergedMFs_idx(alphaSupport, 1);
                else
                    self.mergedIDXs = [ { idxNewMF } self.mergedIDXs ];
                    self.mergedMFs(2:end+1, :) = self.mergedMFs;
                    filtro = compMFs(self.mergedIDXs{2}, self.I_MU) <= newMF(self.I_MU);
                    if any(filtro)
                        self.mergedIDXs{1} = [ self.mergedIDXs{1} ...
                            self.mergedIDXs{2}(filtro) ];
                        self.mergedIDXs{2}(filtro) = [];
                        self.updateMergedMFs_idx(alphaSupport, 2);
                    end
                    self.updateMergedMFs_idx(alphaSupport, 1);
%                     [~, sortMu] = sort(compMFs( self.mergedIDXs{1}, I_MU ));
%                     self.mergedIDXs{1} = self.mergedIDXs{1}(sortMu);
%                     idxSubMF = self.mergedIDXs{2}(1);
%                     while self.possibility( idxSubMF, self.mergedIDXs(2) ) < ...
%                             self.possibility( idxSubMF, self.mergedIDXs(1) )
%                         self.mergedIDXs{2}(1) = [];
%                         self.mergedIDXs{1}(end+1) = idxSubMF;
%                         idxSubMF = self.mergedIDXs{2}(1);
%                     end
                end
            elseif self.mergedMFs(end, self.I_MU) <= newMF(self.I_MU) % se estiver após o último MF
                if self.similarIdx( size(self.mergedMFs, 1), idxNewMF, compMFs ) >= self.Smerge
                    self.mergedIDXs{end}(end+1) = idxNewMF;
                    self.updateMergedMFs_idx(alphaSupport, length(self.mergedIDXs));
                else
                    self.mergedIDXs = [ self.mergedIDXs { idxNewMF } ];
                    self.mergedMFs(end+1, :) = NaN;
                    filtro = newMF(self.I_MU) <= compMFs(self.mergedIDXs{end-1}, self.I_MU);
                    if any(filtro)
                        self.mergedIDXs{end} = [ self.mergedIDXs{end} ...
                            self.mergedIDXs{end-1}(filtro) ];
                        self.mergedIDXs{end-1}(filtro) = [];
                        self.updateMergedMFs_idx(alphaSupport, length(self.mergedIDXs)-1);
                    end
                    self.updateMergedMFs_idx(alphaSupport, length(self.mergedIDXs));
%                     [~, sortMu] = sort(compMFs( self.mergedIDXs{end}, I_MU ));
%                     self.mergedIDXs{end} = self.mergedIDXs{end}(sortMu);
%                     idxSubMF = self.mergedIDXs{end-1}(end);
%                     while self.possibility( idxSubMF, self.mergedIDXs(end-1) ) < ...
%                             self.possibility( idxSubMF, self.mergedIDXs(end) )
%                         self.mergedIDXs{end-1}(end) = [];
%                         self.mergedIDXs{end} = [idxSubMF self.mergedIDXs{end}];
%                         idxSubMF = self.mergedIDXs{end-1}(end);
%                     end
                end
            else
                for i_mf = 2:length(self.mergedIDXs) % para cada MF da feature
                    if  newMF(self.I_MU) < self.mergedMFs(i_mf, self.I_MU) % se estiver antes do MF
                        pre_sim = self.similarIdx( i_mf-1, idxNewMF, compMFs );
                        pos_sim = self.similarIdx( i_mf  , idxNewMF, compMFs );
                        if any([pre_sim, pos_sim] >= self.Smerge)
                            if pre_sim > pos_sim
                                self.mergedIDXs{i_mf-1}(end+1) = idxNewMF;
                                self.updateMergedMFs_idx(alphaSupport, i_mf-1);
                            else
                                self.mergedIDXs{i_mf}(end+1) = idxNewMF;
                                self.updateMergedMFs_idx(alphaSupport, i_mf);
                            end
                        else
                            self.mergedIDXs = [ self.mergedIDXs{1:i_mf-1} ...
                                { idxNewMF }      self.mergedIDXs{i_mf:end} ];
                            self.mergedMFs(i_mf+1:end+1, :) = self.mergedMFs(i_mf:end, :);
                            filtro = newMF(self.I_MU) <= compMFs(self.mergedIDXs{i_mf-1}, self.I_MU);
                            if any(filtro)
                                self.mergedIDXs{i_mf} = [ self.mergedIDXs{i_mf} ...
                                    self.mergedIDXs{i_mf-1}(filtro) ];
                                self.mergedIDXs{i_mf-1}(filtro) = [];
                                self.updateMergedMFs_idx(alphaSupport, i_mf-1);
                            end
                            filtro = compMFs(self.mergedIDXs{i_mf+1}, self.I_MU) <= newMF(self.I_MU);
                            if any(filtro)
                                self.mergedIDXs{i_mf} = [ self.mergedIDXs{i_mf} ...
                                    self.mergedIDXs{i_mf+1}(filtro) ];
                                self.mergedIDXs{i_mf+1}(filtro) = [];
                                self.updateMergedMFs_idx(alphaSupport, i_mf+1);
                            end
                            self.updateMergedMFs_idx(alphaSupport, i_mf);
                            % % %
%                             [~, sortMu] = sort(compMFs( self.mergedIDXs{i_mf-1}, I_MU ));
%                             self.mergedIDXs{i_mf-1} = self.mergedIDXs{i_mf-1}(sortMu);
%                             pre_subMF = self.mergedIDXs{i_mf-1}(end);
%                             [~, sortMu] = sort(compMFs( self.mergedIDXs{i_mf+1}, I_MU ));
%                             self.mergedIDXs{i_mf+1} = self.mergedIDXs{i_mf+1}(sortMu);
%                             pos_subMF = self.mergedIDXs{i_mf+1}(1);
%                             pre_v = self.possibility( pre_subMF, newMF ) ...
%                                 / self.possibility( pre_subMF, self.mergedIDXs(i_mf-1) );
%                             pos_v = self.possibility( pos_subMF, newMF ) ...
%                                 / self.possibility( pos_subMF, self.mergedIDXs(i_mf+1) );
%                             while any([pre_v, pos_v] > 1)
%                                 if pre_v > pos_v
%                                     self.mergedIDXs{i_mf-1} = setdiff(self.mergedIDXs{i_mf-1}, pre_subMF);
%                                     self.mergedIDXs{i_mf}(end+1) = pre_subMF;
%                                     self.mergedMFs(i_mf-1, :) = self.mergeAlphaSupport( alphaSupport( self.mergedIDXs{i_mf-1}, : ) );
%                                     self.mergedMFs(i_mf, :) = self.mergeAlphaSupport( alphaSupport( self.mergedIDXs{i_mf}, : ) );
%                                     pre_subMF = self.mergedIDXs{i_mf-1}(end);
%                                 else
%                                     self.mergedIDXs{i_mf+1} = setdiff(self.mergedIDXs{i_mf+1}, pos_subMF);
%                                     self.mergedIDXs{i_mf}(end+1) = pos_subMF;
%                                     self.mergedMFs(i_mf+1, :) = self.mergeAlphaSupport( alphaSupport( self.mergedIDXs{i_mf+1}, : ) );
%                                     self.mergedMFs(i_mf, :) = self.mergeAlphaSupport( alphaSupport( self.mergedIDXs{i_mf}, : ) );
%                                     pos_subMF = self.mergedIDXs{i_mf+1}(1);
%                                 end
%                                 pre_v = self.possibility( pre_subMF, newMF ) ...
%                                     / self.possibility( pre_subMF, self.mergedIDXs(i_mf-1) );
%                                 pos_v = self.possibility( pos_subMF, newMF ) ...
%                                     / self.possibility( pos_subMF, self.mergedIDXs(i_mf+1) );
%                             end
                        end
                        break;
                    end
                end
            end
        end
        
        function renewMFs(self, compMFs)
            self.mergedIDXs = num2cell(1:size(compMFs, 1));
            alphaSupport = self.updateMergedMFs(compMFs);
            similarities = Inf(length(self.mergedIDXs)-1, 1);
            for ii = 1:length(similarities) % para cada par de componentes vizinhos
                similarities(ii) = self.similarIdx([ii, ii+1], [], compMFs);
            end
            [sim, idx] = max(similarities);
            while sim >= self.Smerge
                self.mergedIDXs{idx} = [self.mergedIDXs{idx} self.mergedIDXs{idx+1}];
                self.updateMergedMFs_idx(alphaSupport, idx);
                self.mergedIDXs(idx+1) = [];
                self.mergedMFs(idx+1, :) = [];
                similarities(idx) = [];
                if idx > 1
                    similarities(idx-1) = self.similarIdx([idx-1, idx], [], compMFs);
                end
                if idx < length(self.mergedIDXs)
                    similarities(idx) = self.similarIdx([idx, idx+1], [], compMFs);
                end
                [sim, idx] = max(similarities);
            end
        end
        
        function purge(self, idx)
            for i_mf = length(self.mergedIDXs):-1:1
                self.mergedIDXs{i_mf} = setdiff(self.mergedIDXs{i_mf}, idx);
                if isempty(self.mergedIDXs{i_mf})
                    self.mergedIDXs(i_mf) = [];
                else
                    filtro = self.mergedIDXs{i_mf} > idx;
                    self.mergedIDXs{i_mf}(filtro) = ...
                        self.mergedIDXs{i_mf}(filtro) - 1;
                end
            end
        end
        
        function [sim, v] = similarIdx(self, idxMergeds, idxMFs, compMFs)
            idxs = [ self.mergedIDXs{idxMergeds} idxMFs ];
            merged = INFGMN_FisVar.mergeAlphaSupport( ...
                INFGMN_FisVar.calcAlphaSupport( ...
                [ self.mergedMFs(idxMergeds, :); compMFs(idxMFs, :) ] ));
            sqrtneglnv = abs( merged(INFGMN_FisVar.I_MU)  - compMFs(idxs, INFGMN_FisVar.I_MU) ) ...
                ./          ( merged(INFGMN_FisVar.I_SPD) + compMFs(idxs, INFGMN_FisVar.I_SPD) );
            [~, idxMinV] = max(sqrtneglnv);
            [sim, v] = self.similarity( merged, compMFs(idxs(idxMinV), :) );
        end
        
        function alphaSupport = updateMergedMFs(self, compMFs)
            alphaSupport = self.calcAlphaSupport(compMFs);
            self.mergedMFs = zeros(length(self.mergedIDXs), 2);
            for ii = 1:length(self.mergedIDXs)
                self.updateMergedMFs_idx(alphaSupport, ii);
            end
        end
        
        function updateMergedMFs_idx(self, alphaSupport, idx)
            self.mergedMFs(idx, :) = self.mergeAlphaSupport( alphaSupport( self.mergedIDXs{idx}, : ) );
        end
        
    end
    
    methods(Static)
        
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
    methods
        
        function [sim, v] = similarity(self, A, B)
            v = self.possibility(A, B);
            if v > (1+self.vmax)/2
                sim = v;
            elseif v < 0.6 % sim < 0.2;
                sim = 0;
            else
                spdmin = min(A(INFGMN_FisVar.I_SPD), B(INFGMN_FisVar.I_SPD));
                beta = 2 * spdmin / (A(INFGMN_FisVar.I_SPD) + B(INFGMN_FisVar.I_SPD));
                sqrtlnv = sqrt(-log(v));
                psi = beta ...
                    + (1 - beta) * erf( (1/(1-beta)) * sqrtlnv ) ...
                    - erf( sqrtlnv );
                sim = psi / (2 - psi);
            end
        end
        
    end
end

