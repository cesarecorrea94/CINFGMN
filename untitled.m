% maxNCs = [20 35 50];
% ndeltas = 5;
% tab = table('Size', [length(maxNCs)*ndeltas 3], ...
%     'VariableTypes', ["double", "double", "double"], ...
%     'VariableNames', {'maxNC',  'delta',  'tau'});
% dim = 2;
% index = 0;
% for deltai = 0:ndeltas
% for maxNC = maxNCs
%     index = index+1;
%     deltamin = 10^(-2-deltai/ndeltas);
%     % dado uma cov mínima (nenhuma componente terá cov menor)
% 	% com NC igualmente espaçadas
%     % no melhor dos casos:
% 	tau = taufordist(...
% 		ones(1,dim)./(2*maxNC^(1/dim)), ... % no ponto menos denso
% 		deltamin, 1)
%     tab(index,:) = {maxNC, deltamin, tau};
% end
% end
% tab


% gmm = INFGMN(minmaxDS(DS), 'delta', 0.0919111184034876, ...
%     'tau', 0.0361844409534152, ...
%     'tmax', 2340, 'spmin', 46.8, 'normalize', false );
% for i=1:10
%     gmm.train(DS);
%     gmm.setFisType('mamdani');
%     y = gmm.recall(DS(:,1));
%     err = dataset2mat(DS(:, end)) - dataset2mat(y(:, end));
%     RMSEml = sqrt(mean(err.^2));
%     gmm.setFisType('sugeno');
%     y = gmm.recall(DS(:,1));
%     err = dataset2mat(DS(:, end)) - dataset2mat(y(:, end));
%     RMSEts = sqrt(mean(err.^2));
%     disp([RMSEml RMSEts]);
% end
% fuzzyLogicDesigner(gmm.toFIS());

%%

% fis=dump_cur.stats(2400).fis;
% f = @() taimesort(fis);
% 
% [   timeit(@tsekouras), timeit(@hefny), timeit(@efsm), timeit(f), timeit(@taimerank); ...
%     timeit(@tsekouras), timeit(@hefny), timeit(@efsm), timeit(f), timeit(@taimerank); ...
%     timeit(@tsekouras), timeit(@hefny), timeit(@efsm), timeit(f), timeit(@taimerank); ...
%     timeit(@tsekouras), timeit(@hefny), timeit(@efsm), timeit(f), timeit(@taimerank) ...
%     ]

% err = zeros(1e1,1e1);
% for ii=1:1e1
%     for jj=ii:1e1
%         A=[ii, 0];
%         B=[jj, 1e1];
%         err(ii,jj) = aproximado(A, B) - similarity(A, B);
%     end
% end
% err

%%

function index = taimesort(fis)
    for inp=1:1e2
        auxtrash = vertcat(fis.Inputs(1).MembershipFunctions.Parameters);
        [~,index] = sortrows(auxtrash(:,2));
        clear auxtrash;
    end
end

function index = taimerank()
    for inp=1:1e2
        index = rankify(randperm(1e2));
    end
end

function rank = rankify(array)
    [sorted, index] = sort(array);
    rank = ones(1, length(array));
    for ii = 2:length(sorted)
        if sorted(ii) == sorted(ii-1)
            rank(ii) = rank(ii-1);
        else, rank(ii) = rank(ii-1)+1;
        end
    end
    rank(index) = rank;
end

%% % % % % % % % % % % % % % %

function tsekouras()
    for inp=1:1e+2
    for ii=1:1e+2
        possibility([ii+1/3 -10], [ii 10]);
        possibility([ii+1/3 -10], [ii 10]);
    end
    end
end
function hefny()
    for inp=1:1e+2
    for ii=1:1e+2
        similarity([ii+1/3 -10], [ii-1/3 10]);
    end
    end
end
function efsm()
    for inp=1:1e+2
    for ii=1:1e+2
        aproximado([ii+1/3 -10], [ii-1/3 10]);
    end
    end
end

function [sim, v] = aproximado(A, B)
    v = possibility(A, B);
    sumAB = (A(1) + B(1))*sqrt(2) * 2*sqrt(-log(0.5));
    base = (A(2) - B(2)) + sumAB;% sumAB = (baseA+baseB)*1 /2 = sumbaseAB /2;
    if base > 0
        % x = [ (x2-x1) *y + (y2x1-y1x2) ]/(y2-y1)
        % y = [(y2-y1)*(y4*x3-y3*x4) - (y4-y3)(y2*x1-y1*x2)] / [(x2-x1)*(y4-y3) - (y2-y1)*(x4-x3)]
        % y = [(A(mu)+A(base)/2) - (B(mu)-B(base)/2)] / [A(base)/2 + B(base)/2)]
        height = base / sumAB;
        intersect = base * height / 2;
        union = sumAB - intersect;
        sim = intersect / union;
    else, sim = 0;
    end
end

function v = possibility(A, B)
    v = exp( -( (A(2)-B(2)) / (A(1)+B(1)) )^2 /2 );
end

function [sim, v] = similarity(A, B)
    spdmin = min(A(1), B(1));
    betasim = 2 * spdmin / (A(1) + B(1));
    v = possibility(A, B);
%     if v > (1+0.8)/2
%         sim = 1;
%     elseif v < 0.6 % sim < 0.2;
%         sim = 0;
%     else
        sqrtlnv = sqrt(-log(v));
        psi = betasim ...
            + (1 - betasim) * erf( (1/(1-betasim)) * sqrtlnv ) ...
            - erf( sqrtlnv );
        sim = psi / (2 - psi);
%     end
end
