classdef INFGMN < handle
    
    properties
        %% Configuration params
        sigma;
        beta;
        delta;
        tau;
        spmin;
        tmax;
        uniform;
        normalize;
        regValue;
        
        sps = [];
        st = [];
        sp = [];
    end
    
    properties (Access = private)
        %% Internal params
        name = mfilename;
        maxDist;
        minDist;
        needsFisUpdate = true;
        
        % Normalization
        proportion;
        minValue = -1;
        maxValue = 1;
        
        % Mahalanobis distance
        mahalaD;
        
        % Components params
        ranges = [];
        priors = [];
        means = [];
        covs = [];
        dill = [];
        nc = 0;
        
        % Sample size
        sampleSize = 0;
        
        % Model likelihood
        dataLikelihood = 0;

        % Components outputs
        loglike = [];
        post = [];
        
        % Fuzzy inference system
        fis;
        inMfType;
        outMfType;
        spread;
        varNames;
        
    end
    
    properties(Constant)
        fisTypes = {'mamdani', 'sugeno'};
        andMethods = {'min', 'prod'};
        orMethods = {'max', 'probor'};
        impMethods = {'min', 'prod'};
        aggMethods = {'max', 'sum', 'probor'};
        defuzzMethods = {'centroid', 'bisector', 'mom', 'lom', 'som', 'wtaver', 'wtsum'};
        mfTypes = {'gaussmf', 'trimf', 'gbellmf', 'trapmf', 'pimf', 'psigmf', 'dsigmf'};
        defaultFISOptions = struct(...
            'name', mfilename, ...
            'type', 'mamdani', ...
            'andMethod', 'prod', ...
            'orMethod', 'max', ...
            'impMethod', 'prod', ...
            'aggMethod', 'sum', ...
            'defuzzMethod', 'centroid',...
            'inMfType', 'gaussmf', ...
            'outMfType', 'gaussmf', ...
            'spread', 0.5);
    end
    
    methods
        
        %% Constructor
        function self = INFGMN(ranges, varargin)
            
            
            validateattributes(ranges, {'dataset'}, {'nonempty', 'nrows', 2}, self.name, '"ranges"')
            
            self.ranges = double(ranges);
            self.varNames = ranges.Properties.VarNames;
            dimension = size(self.ranges, 2);
            
            parser = inputParser;
            parser.StructExpand = false;
            
            defaultDelta = 0.01;
            defaultBeta = 1;
            defaultTau = 0.1;
            defaultTMax = 2 * dimension + 1;
            defaultSPMin = dimension + 1;
            defaultUniform = false;
            defaultNormalize = true;
            defaultRegValue = 0;
            
            checkPositiveInterval = @(x, pname) validateattributes(x, {'numeric'}, {'>=', 0,'<=',1, 'scalar'}, self.name, pname);
            checkSigma = @(x) validateattributes(x, {'numeric'}, {'nonempty', 'row', '>', 0}, self.name, '"sigma"');
            checkLogical = @(x, pname) validateattributes(x, {'logical'}, {'nonempty'}, self.name, pname);
            checkPositiveNumScalar = @(x, pname) validateattributes(x, {'numeric'}, {'>', 0, 'scalar'}, self.name, pname);
            checkFISOptions = @(x) validateattributes(x, {'struct'}, {'nonempty'}, self.name, '"FIS_Options"');
            
            addOptional(parser, 'sigma', [], checkSigma);
            addOptional(parser, 'delta', defaultDelta, @(x) checkPositiveInterval(x, '"delta"'));
            addOptional(parser, 'beta', defaultBeta, @(x) checkPositiveInterval(x, '"beta"'));
            addOptional(parser, 'tau', defaultTau, @(x) checkPositiveInterval(x, '"tau"'));
            addOptional(parser, 'tmax', defaultTMax, @(x)checkPositiveNumScalar(x, '"tmax"'));
            addOptional(parser, 'spmin', defaultSPMin, @(x)checkPositiveNumScalar(x, '"spmin"'));
            addOptional(parser, 'regvalue', defaultRegValue, @(x)checkPositiveNumScalar(x, '"regvalue"'));
            addOptional(parser, 'uniform', defaultUniform, @(x) checkLogical(x, 'uniform'));
            addOptional(parser, 'normalize', defaultNormalize, @(x) checkLogical(x, 'normalize'));
            addOptional(parser, 'FIS_Options', self.defaultFISOptions, checkFISOptions);
            
            parse(parser, varargin{:});
            
            self.normalize = parser.Results.normalize;
            self.sigma = parser.Results.sigma;
            
            if self.normalize
                [ranges, self.proportion] = mapminmax(self.ranges', self.minValue, self.maxValue);
                self.ranges = ranges';
            end
            
            rangeLengths = abs(max(self.ranges) - min(self.ranges));
            self.delta = parser.Results.delta;
            
            self.dill = rangeLengths * 0.001;
            
            if isempty(self.sigma) && (self.delta ~= defaultDelta)
                self.sigma = (self.delta * rangeLengths) .^ 2;
            else
                self.sigma = (defaultDelta * rangeLengths) .^ 2;
            end
            
            if size(self.sigma, 2) ~= length(rangeLengths)
                throw(MException(['MATLAB:' self.name ':expectedEqualColumnNumber'], ...
                    '"ranges" and "sigma" must have the same number of columns'));
            end
            
            self.tmax = parser.Results.tmax;
            self.beta = parser.Results.beta;
            self.spmin = parser.Results.spmin;
            self.uniform = parser.Results.uniform;
            self.regValue = parser.Results.regvalue;
            self.tau = parser.Results.tau;
            self.maxDist = chi2inv(1 - self.tau, dimension);
            
            userOptions = parser.Results.FIS_Options;
            fn = fieldnames(userOptions);
            options = self.defaultFISOptions;
            for i = 1:length(fn)              
                options.(fn{i}) = userOptions.(fn{i});
            end
            
            optionsParser = inputParser;

            checkFISName = @(x) validateattributes(x, {'char'}, {'nonempty'}, self.name, '"FIS_Options.name"');
            checkFISType = @(x) any(validatestring(x, self.fisTypes));
            checkFISMethods = @(x, validMethods, optionName) ... 
                (ischar(x) && exist(x) == 2) || any(validatestring(x, validMethods, self.name, optionName)); %#ok<*EXIST>
            
            checkMfType = @(x, pname) any(validatestring(x, self.mfTypes, self.name, pname)) ;
            
            checkSpread = @(x, pname) validateattributes(x, {'numeric'}, {'scalar', '>', 0, '<=', 1});

            addParameter(optionsParser,'andMethod', [], @(x) checkFISMethods(x, self.andMethods, '"FIS_Options.andMethod"'));
            addParameter(optionsParser,'orMethod', [], @(x) checkFISMethods(x, self.orMethods, '"FIS_Options.orMethod"'));
            addParameter(optionsParser,'impMethod', [], @(x) checkFISMethods(x, self.impMethods, '"FIS_Options.impMethod"'));
            addParameter(optionsParser,'aggMethod', [], @(x) checkFISMethods(x, self.aggMethods, '"FIS_Options.aggMethod"'));
            addParameter(optionsParser,'defuzzMethod', [], @(x) checkFISMethods(x, self.defuzzMethods, '"FIS_Options.defuzzMethod"'));
            addParameter(optionsParser,'name', [], checkFISName);
            addParameter(optionsParser,'type', [], checkFISType);
            addParameter(optionsParser,'inMfType', [], @(x) checkMfType(x, '"FIS_Options.inMfType"'));
            addParameter(optionsParser,'outMfType', [], @(x) checkMfType(x, '"FIS_Options.outMfType"'));
            addParameter(optionsParser,'spread', [], @(x) checkSpread(x, '"FIS_Options.spread"'));
            
            parse(optionsParser, options);
            
            %% Create the fuzzy layer
            self.name = optionsParser.Results.name;
            self.spread = optionsParser.Results.spread;
            type = optionsParser.Results.type;
            andMethod = optionsParser.Results.andMethod;
            orMethod = optionsParser.Results.orMethod;
            impMethod = optionsParser.Results.impMethod;
            aggMethod = optionsParser.Results.aggMethod;
            defuzzMethod = optionsParser.Results.defuzzMethod;
    
            self.fis = newfis(self.name, type, andMethod, orMethod, impMethod, aggMethod, defuzzMethod);
            self.inMfType = optionsParser.Results.inMfType;
            self.outMfType = optionsParser.Results.outMfType;
            
            %% Creates a component at the mean of the estimated ranges of the problem features.
            self.createComponent(mean(self.ranges))
        end
        
         %% Train the INFGMN with the data set X
        function self = train(self, X)
            
            validateattributes(X, {'dataset'}, {'nonempty', 'ncols', size(self.ranges, 2)}, strcat(self.name, ':train'),'"X"');
            X = double(X);
            
            if self.normalize
                X = mapminmax('apply', X', self.proportion)';
            end
            
            if any(any(isnan(X)))
                X = self.internal_imputation(X);
            end
            
            N = size(X,1);
            for i = 1:N
                x = X(i,:);
                self.computeLikelihood(x);
                if (~self.hasAcceptableDistribution())
                    self.createComponent(x);
                end
                self.computePosterior();
                self.updateComponents(x);
                self.sampleSize = self.sampleSize + 1;
                self.removeSpurious();
            end
            
            % Force fuzzy layer update.
            self.needsFisUpdate = true;
%             disp(self.nc);
        end
        

        %% Compute the log probability density of the multivariate normal distribution
        % evaluated for each INFGMN gaussian component.
        function computeLikelihood(self, x)
            if self.nc > 0
                [self.loglike, self.mahalaD] = self.logmvnpdf(x, self.means, self.covs);
            end
        end
        
        %% Compute the log probability density of a multivariate normal distribution
        function [loglike, mahalaD] = logmvnpdf(self, x, means, covs)
            [n,d] = size(x);
            k = size(means,1);
            loglike = zeros(n, k);
            mahalaD = zeros(n, k);
            logSqrtDetCov = zeros(1, k);

            for j = 1:k
                [L,err] = cholcov(covs(:,:,j), 0);
                if err ~= 0
%                     error('stats:mvnpdf:BadCovariance',...
%                         'Each page of SIGMA must be symmetric and positive definite.');
                    L = sqrt(self.jdiag(covs(:,:,j))); % a vector
                    if  any(L < eps(max(L)) * d)
                        error('Ill Conditioned Covariance.');
                    end
                    L = diag(L); % a matrix
                end
                logSqrtDetCov(j) = sum(log(self.jdiag(L)));

                Xcentered = bsxfun(@minus, x, means(j,:));
%                 xRinv = bsxfun(@times, Xcentered , (1 ./ L));
                xRinv = Xcentered / L;

                mahalaD(:,j) = sum(xRinv.^2, 2);
                loglike(:,j) = - 0.5 * (mahalaD(:,j) + 2*logSqrtDetCov(j) + d * log(2 * pi));
            end
        end
  

        %% Compute the posterior probability for each INFGMN gaussian component.
        function computePosterior(self)
            logPrior = log(self.priors);
            self.post = bsxfun(@plus, self.loglike, logPrior);
            maxll = max(self.post, [], 2);
            % minus maxll to avoid underflow
            self.post = exp(bsxfun(@minus, self.post, maxll));
            density = sum(self.post, 2);
            % normalize
            self.post = bsxfun(@rdivide, self.post, density);
            logpdf = log(density) + maxll;
            self.dataLikelihood = self.dataLikelihood + sum(logpdf);
        end

        %% Check the INFGMN novelty criterion.
        function h = hasAcceptableDistribution(self)
            for j = 1:self.nc
                if (self.mahalaD(j) <= self.maxDist)
                    h = true;
                    return;
                end
            end
            h = false;
        end

        %% Create a gaussian component when a data point x matches the novelty criterion.
        function createComponent(self, x)
            self.nc = self.nc + 1;
            self.means(self.nc,:) = x;
            self.sps(1,self.nc) = 1;
            self.st(1,self.nc) = 0;
            self.sp(1,self.nc) = 0;
            self.updatePriors();
            if self.nc > 1
                self.covs(:,:,self.nc) = diag(mean(self.jdiag(self.covs), 3));
            else
                self.covs(:,:,self.nc) = diag(self.sigma);
            end
            [self.loglike(1,self.nc), self.mahalaD(1, self.nc)] = ...
                self.logmvnpdf(x, self.means(self.nc,:), self.covs(:,:,self.nc));
        end

        %% Update the INFGMN gaussians priors.
        function updatePriors(self)
            if ~self.uniform
                self.priors = self.sps ./ sum(self.sps);
             else
                self.priors = ones(1,self.nc) ./ self.nc;
            end
        end

        %% Update the INFGMN components parameters.
        function updateComponents(self, x)
            self.sps = self.sps + self.post;
            self.st = self.st + 1.0;
            self.sp = self.sp + self.post;
            w = (self.post ./ self.sps) * self.beta;
            for j = 1:self.nc
                Xcentered = x - self.means(j,:);
                deltaMU = w(j) * Xcentered;
                self.means(j,:) = self.means(j,:) + deltaMU;
                XcenteredNEW = x - self.means(j,:);
                self.covs(:,:,j) = self.covs(:,:,j) - (deltaMU')*(deltaMU) + ...
                    w(j) * ((XcenteredNEW')*(XcenteredNEW) - self.covs(:,:,j)) + self.regValue; 
                for d = 1:length(self.varNames)
                    if self.covs(d,d,j) < 0
                        self.covs(d,d,j) = 0.0000001;
                    end
                end
            end
%             self.covs(self.covs < 0) = 0.0000001;
            % normalize the priors
            self.updatePriors();
        end
     
        %% Remove spurious gaussian components.
        function removeSpurious(self)
            for i = self.nc:-1:1
                if (self.st(i) == self.tmax)
                    if self.sp(i) < self.spmin
                        self.nc = self.nc - 1;
                        self.priors(i) = [];
                        self.means(i,:) = [];
                        self.covs(:,:,i) = [];
                        self.st(i) = [];
                        self.sps(i) = [];
                        self.sp(i) = [];
                        self.post(i) = [];
                        self.mahalaD(i) = [];
                        self.loglike(i) = [];
                    else
                        self.st(i) = 0;
                        self.sp(i) = 0;
                    end
                end
            end
        end

        %% INFGMN recalling algorithm
        function y = recall(self, x)
            validateattributes(x, {'dataset'}, {'nonempty', }, strcat(self.name, ':recall'),'"x"');
            varNames = x.Properties.VarNames; %#ok<PROPLC>
            inputIndexes = cellfun( @(var)find(strcmp(self.varNames, var)), varNames, 'UniformOutput', 1); %#ok<PROPLC>
            outputIndexes = 1:length(self.varNames);
            outputIndexes(inputIndexes) = []; 
            self.updateFuzzyLayer(inputIndexes, outputIndexes);
            if self.normalize
                y = ones(size(x, 1), length(self.varNames));
                y(:, inputIndexes) = double(x);
                y = mapminmax('apply', y', self.proportion)';
                if any(any(isnan(y(:, inputIndexes))))
                    x_imp = nan(size(x, 1), length(self.varNames));
                    x_imp(:, inputIndexes) = y(:, inputIndexes);
                    x_imp = self.internal_imputation(x_imp);
                    y(:, inputIndexes) = x_imp(:, inputIndexes); 
                end
                y(:, outputIndexes) = evalfis(y(:, inputIndexes), self.fis);
                y = mapminmax('reverse', y', self.proportion)';
                y = y(:, outputIndexes);
            else
                y = double(x);
                if any(any(isnan(y)))
                    x_imp = nan(size(x, 1), length(self.varNames));
                    x_imp(:, inputIndexes) = y;
                    x_imp = self.internal_imputation(x_imp);
                    y = x_imp(:, inputIndexes);
                end
                y = evalfis(y, self.fis);
            end
            y = mat2dataset(y , 'VarNames', self.varNames(outputIndexes));
        end
            
        % INFGMN internal_imputation method
        function x = internal_imputation(self, x)
            indexes = isnan(x);
            if all(indexes)
                error('Each trainning input pattern must have at least 1 non nan column.');
            end;  
            if any(any(indexes))
                N = size(x, 1);
                for j = 1:N
                    pajs = zeros(self.nc, 1);
                    curIndexes = indexes(j, :);
                    xm = zeros(self.nc, sum(curIndexes));

                    for i = 1:self.nc
                        meanA = self.means(i, ~curIndexes);
                        meanB = self.means(i, curIndexes);
                        covA = self.covs(~curIndexes, ~curIndexes, i);
                        xm(i,:) = meanB;
                        ll = self.logmvnpdf(x(j, ~curIndexes), meanA, covA);
                        pajs(i) = exp(ll) * self.priors(i);
                    end
                    pajs = pajs ./ sum(pajs);
                    x(j, curIndexes) = sum(bsxfun(@times, xm, pajs));
                    x(j, x(j, :) > self.ranges(2, :)) = self.ranges(2, x(j, :) > self.ranges(2, :));
                    x(j, x(j, :) < self.ranges(1, :)) = self.ranges(1, x(j, :) < self.ranges(1, :));
                end
            end        
        end
        
        % INFGMN imputation method
        function y = imputation(self, x)
            validateattributes(x, {'dataset'}, {'nonempty', 'nrows', 1}, strcat(self.name, ':imputation'),'"x"');
            values = double(x);
            indexes = isnan(values);
            if all(indexes)
                error('Each trainning input pattern must have at least 1 non nan column.');
            end;
            if any(any(indexes))
                y = ones(1, length(self.varNames));
                y(~indexes) = values(~indexes);
                if self.normalize
                    y = mapminmax('apply', y', self.proportion)';
                end
                pajs = zeros(self.nc, 1);
                xm = zeros(self.nc, sum(indexes));

                for i = 1:self.nc
                    meanA = self.means(i, ~indexes);
                    meanB = self.means(i, indexes);
                    covA = self.covs(~indexes, ~indexes, i);
                    xm(i,:) = meanB;
                    ll = self.logmvnpdf(y(~indexes), meanA, covA);
                    pajs(i) = exp(ll) * self.priors(i);
                end
                pajs = pajs ./ sum(pajs);
                y(indexes) = sum(bsxfun(@times, xm, pajs));
                if self.normalize
                    y = mapminmax('reverse', y', self.proportion)';             
                end
                y = mat2dataset(y(indexes) , 'VarNames', self.varNames(indexes));
            end
        end
        
        % Calculate the indexes assignments for the data points in X for
        % cluster analysis.
        function idx = cluster(self, X)
            N = size(X,1);
            idx = zeros(N,1);
            for i=1:N
                x = X(i,:);
                self.computeLikelihood(x);
                self.computePosterior();
                [~,idx(i)] = max(self.post);
            end
        end
        
        
        %% Accessor methods.
        function ms = modelSize(self)
            %ms = size(self.fis.rule, 2);
            ms = self.nc;
        end
        
        function covs = modelCovariances(self)
            covs = self.covs;
        end
        
        function means = modelMeans(self)
            means = self.means;
        end
        
        function weights = modelWeights(self)
            weights = self.priors;
        end
        
        function spreads = modelSpreads(self)
            if size(self.covs, 3) > 1
                spreads = squeeze(self.jdiag(self.covs))' .^ self.spread;
            else
                spreads = squeeze(self.jdiag(self.covs)) .^ self.spread;
            end
        end
        
        function fis = toFIS(self)
            fis = self.fis;
        end
        
        function mfTypes = getMfTypes(self)
            mfTypes = {self.inMfType, self.outMfType};
        end

        function diagcovs = jdiag(~, covs)
            diagcovs = zeros(1, size(covs,2), size(covs,3));
            for j = 1:size(covs,3)
                diagcovs(:,:,j) = diag(covs(:,:,j))';
            end
        end
        
        function setFisType(self, newFisType)
            if strcmp(newFisType, 'mamdani')
                self.outMfType = 'gaussmf';
                defuzzMethod = 'centroid';
            elseif strcmp(newFisType, 'sugeno')
                self.outMfType = 'linear';
                defuzzMethod = 'wtaver';
            else
                throw(MException(['MATLAB:' self.name ':incorrectFisType'], ...
                    '"newFisType" must be "mamdani" or "sugeno"'));
            end
            self.fis = newfis(self.name, newFisType, 'prod', 'max', 'prod', 'sum', defuzzMethod);
            self.needsFisUpdate = true;
        end
        
        function setSpread(self, newSpread)
            self.spread = newSpread;
            self.needsFisUpdate = true;
        end
        
    end %end methods
    
    methods (Access = private)
        
        function updateFuzzyLayer(self, inputIndexes, outputIndexes)
            self.needsFisUpdate = self.needsFisUpdate || ...
                isempty(self.fis.output) || ...
                ~isempty(setdiff([self.fis.output(:).name], self.varNames(outputIndexes))) || ...
                ~isempty(setdiff(self.varNames(outputIndexes),[self.fis.output(:).name]));
            if (self.needsFisUpdate)
                self.cleanFis();
                self.createFuzzyVariables(inputIndexes, outputIndexes);
                self.createRules();
            end
            self.needsFisUpdate = false;
        end
        
        function createFuzzyVariables(self, inputIndexes, outputIndexes)
            for i = 1:length(inputIndexes)
                range = self.ranges(:, inputIndexes(i));
%                 mus = self.dilatedMeansForVar(inputIndexes(i)); 
                mus = self.means(:, inputIndexes(i));
                spreads = self.jdiag(self.covs(inputIndexes(i),inputIndexes(i),:)) .^ self.spread;
                self.fis = addvar(self.fis, 'input', ...
                    self.varNames{inputIndexes(i)}, [range(1), range(2)]);
                for j = 1:self.nc
                    mu = mus(j);
                    spd = spreads(j);
                    params = self.toMfParams('trimf', mu, spd, self.inMfType);
                    if self.fis.input(i).range(1) >  params(1)
                        self.fis.input(i).range(1) = params(1) - 0.3;
                    else
                        if self.fis.input(i).range(2) <  params(3)
                            self.fis.input(i).range(2) = params(3) + 0.3;
                        end
                    end
                    self.addInputMf(i, mu, spd, ['MF' num2str(j)])
                end
            end
            for i = 1:length(outputIndexes)
                range = self.ranges(:, outputIndexes(i));
                self.fis = addvar(self.fis, 'output', ...
                    self.varNames{outputIndexes(i)}, [range(1), range(2)]);
            end
            if strcmp(self.fis.type, 'mamdani')
                for i = 1:length(outputIndexes)
%                     mus = self.dilatedMeansForVar(outputIndexes(i));
                    mus = self.means(:, outputIndexes(i));
                    spreads = self.jdiag(self.covs(outputIndexes(i),outputIndexes(i),:)) .^ self.spread;
                    for j = 1:self.nc
                        mu = mus(j);
                        spd = spreads(j);
                        params = self.toMfParams('trimf', mu, spd, self.outMfType);
                        if self.fis.output(i).range(1) >  params(1)
                            self.fis.output(i).range(1) = params(1) - 0.3;
                        else
                            if self.fis.output(i).range(2) <  params(3)
                                self.fis.output(i).range(2) = params(3) + 0.3;
                            end
                        end
                        self.addOutputMf(i, mu, spd, ['MF' num2str(j)])
                    end
                end
            else % sugeno
                if ~strcmp(self.outMfType, 'linear')
                    throw(MException(['MATLAB:' self.name ':errorOnSugenoOutputType'], ...
                        'outputMF type must be linear (constant type not supported yet)'));
                end
                for j = 1:self.nc
                    muj = self.means(j,:);
                    invcovj = pinv(sign(self.covs(:,:,j)) ...
                        .* abs(self.covs(:,:,j)) .^ self.spread);
                    for i = 1:length(outputIndexes)
                        o = outputIndexes(i);
                        CAngular = -invcovj(inputIndexes, o) / invcovj(o, o);
                        CLinear = muj(o) - muj(inputIndexes) * CAngular;
                        params = [CAngular' CLinear];
                        self.fis = addmf(self.fis, 'output', i, ...
                            ['MF' num2str(j)], self.outMfType, params);
                    end
                end
            end
        end
        
        function createRules(self)
            numVars = length(self.varNames);
            mfsTemplate = ones(1, numVars);
            ruleList = ones(self.nc, numVars + 2);
            for i = 1:self.nc  
                ruleList(i, :) = [(mfsTemplate .* i) self.priors(i) 1];
            end
            self.fis = addrule(self.fis, ruleList);
        end
        
        function cleanFis(self)
            self.fis.rule = [];
            self.fis.output = [];
            self.fis.input = [];
        end
        
        function addInputMf(self, index, mean, spread, mfName)
            params = self.toMfParams(self.inMfType, mean, spread, self.defaultFISOptions.inMfType);
            self.fis = addmf(self.fis, 'input', index, mfName, self.inMfType, params);
        end
        
        function addOutputMf(self, index, mean, spread, mfName)
            params = self.toMfParams(self.outMfType, mean, spread, self.defaultFISOptions.outMfType);
            self.fis = addmf(self.fis, 'output', index, mfName, self.outMfType, params);
        end
        
        function params = toMfParams(~, mfType, mean, spread, defaultMfType)
            params = [spread mean];
            if ~strcmp(mfType, defaultMfType)
                params = mf2mf(params, defaultMfType, mfType);
            end
        end 
        
        function dilMeans = dilatedMeansForVar(self, varIndex)
            dilMeans = self.means(:, varIndex);
            if self.nc > 1
                dil = self.dill(varIndex);
                [sorted, indexes] = sort(dilMeans, 'ascend'); 
                dilMeans(indexes(end)) = sorted(end) + dil;
                dilMeans(indexes(1)) = sorted(1) - dil;
    
            end
        end
        
        function dump(self)
            cvs = self.modelSpreads();
            numCols = size(self.means, 2);            
            disp('================================================');
            disp('rules:')
            fprintf('%7s\n', num2str(self.nc));
            disp('------------------------------------------------');
            disp('means:');
            fprintf('%s', sprintf([repmat('%16.8f', 1, numCols) '\n'], self.means'));
            disp('------------------------------------------------');
            disp('covs:');
            fprintf('%s', sprintf([repmat('%16.8f', 1, numCols) '\n'], cvs'));
            disp('------------------------------------------------');
            disp('st:');
            fprintf('%s', sprintf([repmat('%16.8f', 1, numCols) '\n'], self.st));
            if (mod(length(self.st), 3) ~= 0); disp(' '); end
            disp('------------------------------------------------');
            disp('sp:');
            fprintf('%s', sprintf([repmat('%16.8f', 1, numCols) '\n'], self.sp));
            if (mod(length(self.sp), 3) ~= 0); disp(' '); end
            disp('------------------------------------------------');
            disp('sps:');
            fprintf('%s', sprintf([repmat('%16.8f', 1, numCols) '\n'], self.sps));
            if mod(length(self.sps), 3) ~= 0; disp(' '); end
            disp('------------------------------------------------');
            disp('mahalaD:');
            fprintf('%s', sprintf([repmat('%16.8f', 1, numCols) '\n'], self.mahalaD));
            if mod(length(self.mahalaD), 3) ~= 0; disp(' '); end
            disp('------------------------------------------------');
            disp('loglike:');
            fprintf('%s', sprintf([repmat('%16.8f', 1, numCols) '\n'], self.loglike));
            if mod(length(self.loglike), 3) ~= 0; disp(' '); end
            disp('------------------------------------------------');
            disp('priors:');
            fprintf('%s', sprintf([repmat('%16.8f', 1, numCols) '\n'], self.priors));
            if mod(length(self.priors), 3) ~= 0; disp(' '); end
            disp('------------------------------------------------');
            disp('post:');
            fprintf('%s', sprintf([repmat('%16.8f', 1, numCols) '\n'], self.post));
            if mod(length(self.post), 3) ~= 0; disp(' '); end
            disp('================================================');
        end
    end
    
end % end class


