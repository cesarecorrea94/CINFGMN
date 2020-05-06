
function gmm = sinsample(delta, tau, tmax, spmin, normalize)
    warning off;
    rng(0);
    DS = sin_DS(true);
    gmm = INFGMN(minmaxDS(DS), 'delta', delta, 'tau', tau, ...
        'tmax', tmax, 'spmin', spmin, 'normalize', normalize );
    gmm.train(DS);
%     fuzzyLogicDesigner(gmm.toFIS());

    DS = sin_DS(false);
    actual = dataset2mat(DS(:, end));
    % 
    gmm.setFisType('mamdani');
    y = gmm.recall(DS(:,1));
    mamdani_output = dataset2mat(y(:, end));
    mamdani_err = abs(actual - mamdani_output);
    % 
    gmm.setFisType('sugeno');
    y = gmm.recall(DS(:,1));
    sugeno_output = dataset2mat(y(:, end));
    sugeno_err = abs(actual - sugeno_output);
    % 
    fig = figure;
    tiledlayout(2,1)
    % % % % % 
    ax1 = nexttile;%newaxis(fig);
    hold(ax1, 'on');

    plot(DS(:,1),       actual,    'k-', 'DisplayName', 'actual');
    plot(DS(:,1), mamdani_output, 'b-.', 'DisplayName', 'ML output', 'LineWidth', 1.2);
    plot(DS(:,1),  sugeno_output, 'r--', 'DisplayName', 'TS output', 'LineWidth', 1.2);

    plotlabels(ax1, 'x', 'sin(x)');
    legend;
    % % % % % 
    ax1 = nexttile;%newaxis(fig);
    hold(ax1, 'on');

    plot(DS(:,1), mamdani_err, 'b-.', 'DisplayName', 'ML error', 'LineWidth', 1.2);
    plot(DS(:,1),  sugeno_err, 'r--', 'DisplayName', 'TS error', 'LineWidth', 1.2);

    plotlabels(ax1, 'x', 'abs(error)');
    legend;
    % 
    saveas(fig, ['imgs/sin', ...
        '_normalize', num2str( normalize ), ...
        '_rules', num2str( gmm.modelSize() ) ...
        ], 'png');
    close(fig);
end
