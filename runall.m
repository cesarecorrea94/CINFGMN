warning off;

% try
%     abbaalone
% catch ME
%     disp(ME);
% end
% clear;
% 
% try
%     concreto
% catch ME
%     disp(ME);
% end
% clear;
% 
% try
%     housing
% catch ME
%     disp(ME);
% end
% clear;
% 
% try
%     slump
% catch ME
%     disp(ME);
% end
% clear;

%%

uci = { 'abbaalone', 'concreto', 'housing', 'slump' };
% for idxuci = 1:4
%     dumpname = ['dumps/' uci{idxuci} '/29-Jun-2020 filter.mat'];
%     dumpseries = INFGMN_series(dumpname);
%     self = dumpseries.myself();
%     self.save_fis = true;
%     [self.dumps.cputime] = deal(nan);
%     [self.dumps.stats] = deal([]);
%     assert(self.doMerge);
%     assert(self.recallTest);
%     save([dumpseries.dumpname '.merge_test.mat'], '-struct', 'self');
% 
%     self.recallTest = false;
%     save([dumpseries.dumpname '.merge_train.mat'], '-struct', 'self');
% 
%     self.doMerge = false;
%     self.fis_types = {'mamdani', 'sugeno'};
%     self.recallTest = true;
%     save([dumpseries.dumpname '.dont_test.mat'], '-struct', 'self');
% 
%     self.recallTest = false;
%     save([dumpseries.dumpname '.dont_train.mat'], '-struct', 'self');
% end

ramos = {'merge_test', 'merge_train', ...
          'dont_test',  'dont_train'};
% for idxuci = 1:4
%     dumpname = ['dumps/' uci{idxuci} '/29-Jun-2020 filter.mat'];
%     for idxramos = 1:4
%         dumpseries = INFGMN_series([dumpname '.' ramos{idxramos} '.mat']);
%         dumpseries.update(@INFGMN_series.step_nonseries);
%     end
% end
% 
% error('Finished');

% ramos = {'merge_test', 'dont_test'};
for idxuci = 4:4%1:4
    dumpname = ['dumps/' uci{idxuci} '/29-Jun-2020 filter.mat'];
    for idxramos = 1:length(ramos)
        dumpseries = INFGMN_series([dumpname '.' ramos{idxramos} '.mat']);
        self = dumpseries.myself();
        eval([ramos{idxramos} ' = self.dumps;']);
%         print_all_fis(dumpseries.myself(), 2, ...
%             ['imgs/' uci{idxuci} '/29-Jun ' ramos{idxramos}]);
    end
end
% clear;
