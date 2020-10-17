
function print_all_fis(self, topN, folder)
    ranges = minmaxDS(self.DS);
    for ii = 1:topN
        subfolder = [folder '/' num2str(ii)];
        stats = self.dumps(ii).stats;
        for jj = 1:length(stats)
            if  isnumeric(stats(jj).fis)
                continue;
            end
            mkdir(subfolder);
            myplotfis(stats(jj).fis, 1, ranges, ...
                ...[subfolder '/' num2str(jj) '_NC' num2str(stats(jj).NC)]); % UCI
                [subfolder '/t' num2str(stats(jj).t)]); % DynNLS
        end
    end
end
