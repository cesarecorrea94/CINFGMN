
function getINFGMNParam(self)
    log2spmin = self.dumps.log2tmax - self.dumps.log2maxNC;
    disp({'maxNC',      2^ [self.dumps.log2maxNC]  ;
        'delta',        2^ [self.dumps.log2delta]  ;
        'tau',          2^ [self.dumps.log2tau]    ;
        'tmax', round(  2^ [self.dumps.log2tmax] ) ;
        'spmin',        2^          log2spmin      ;
        'doMerge', self.doMerge ;
        'simMerge', self.Smerge ;
        'simRefit', NaN;%self.Sdeath ;
        'maxMF', self.maxFoCSize });
end
