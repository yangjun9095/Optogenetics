% Plot the NLS-mCherry-LEXY export for different excitation wavelength

% Load the dataset
%File directories
FigPath = 'S:\YangJoon\Dropbox\Optogenetics\Figures';
DataPath = 'S:\YangJoon\Dropbox\Optogenetics\WavelengthDependence';
D = dir(DataPath);

% Note. D is a cell, I need to use cell2mat to extract from D
data_compiled = {};
k=1; % counter

for i=3:length(D)
    % extract the wavelength from the file name
    filename = D(i).name;
    lambda = filename(end-2:end);
    
    % define a data structure to extract compilednuclei info.
    dataset = struct;
    dataset = load([DataPath, filesep, D(i).name, filesep, 'CompiledNuclei.mat']);
    
    % extract key fields for generating plots
    MeanFluo = dataset.MeanVectorAll;
    SDFluo = dataset.SDVectorAll;
    SEFluo = dataset.SDVectorAll./sqrt(dataset.NParticlesAll);
    
    Time = dataset.ElapsedTime;
    Frame = 1:length(MeanFluo);

    data_compiled(k).wavelength = str2num(lambda);
    data_compiled(k).MeanFluo = MeanFluo;
    data_compiled(k).SEFluo = SEFluo;
    data_compiled(k).Time = Time;
    data_compiled(k).Frames = Frame;
    
    k=k+1; % counter count
end


%% generate plots of time traces

colorCode = [];
wavelength = {};
k=1; % counter

hold on
for i= [1,3,5,8,10]%1:length(data_compiled)
    % Normalization of the fluorescence intensity is needed 
    % as the laser power is not controlled in this case.
    % We tried with strong enough laser power, to saturate the export.
    
    % color code matching the wavelength
    colorCode(i,:) = wavelength2color(data_compiled(i).wavelength + 40*k, 'colorSpace', 'rgb');
    
    % errorbar plot of fluo intensity over time
    MaxIntensity = max(data_compiled(i).MeanFluo);
    % find the time point where it's the maximum
    MaxIndex = find(data_compiled(i).MeanFluo == MaxIntensity);
    
    errorbar(data_compiled(i).Time,...
                data_compiled(i).MeanFluo/MaxIntensity,...
                data_compiled(i).SEFluo/MaxIntensity)%,...
                %'Color',colorCode(i,:))
    wavelength{k} = num2str(data_compiled(i).wavelength);
    k=k+1;
end

xlim([0 5])
ylim([0 1.2])
xlabel('time (min)')
ylabel('normalized intensity')
xticks([0 1 2 3 4 5])
yticks([0 0.5 1 1.5])
legend(wavelength)
box on
StandardFigure(gcf,gca)

saveas(gcf, [FigPath, filesep, 'NLS-mCherry-LEXY-lambda_dependence_Maxnormalized.pdf'])
%% Test the color coding function
wave = 458:1:670;
hold on
for i=1:length(wave)
    colorCode = wavelength2color(wave(i), 'colorSpace', 'rgb')
    plot(1:10, (1:10)*i/100,'Color',colorCode)
end

%% plot to check (time window for the fitting)
hold on
for i=1:length(data_compiled)
    Time = [];
    Fluo = [];
    Fluo_SE = [];
    
    Time = data_compiled(i).Time;
    Fluo = data_compiled(i).MeanFluo;
    Fluo_SE = data_compiled(i).SEFluo;
    errorbar(Time, Fluo/(max(Fluo)), Fluo_SE/(max(Fluo)))
%     plot(Time, log(Fluo/(max(Fluo))))
    %pause
end
%% fitting an exponential decay function to extract the decay constant
for i=1:length(data_compiled)
    % extract the useful fields, time and fluo
    Time = [];
    Fluo = [];
    
    Time = data_compiled(i).Time;
    Fluo = data_compiled(i).MeanFluo;

    % find the time point where it's the maximum
    MaxFluo = max(Fluo(5:end)); % this is because we illuminated from the end of the 5th frame.
    if i==10
        MaxIndex = 7;
    elseif i==3 || i==4 || i==5
        MaxIndex = 5;
    else
        MaxIndex = 6;%find(Fluo == MaxFluo);
    end
    
    % normalize the fluo to its maximum
    Fluo = Fluo/MaxFluo;
    
    % truncate the data from the maximum fluorescence to focus on the decay
    % regime.
    Time = Time(MaxIndex:end) - Time(MaxIndex);
    Fluo = Fluo(MaxIndex:end);
    decay_model = @(params) params(1)*exp(-params(2)*Time) + params(3) - Fluo; 
    
    % initial guess for the parameters
    params0 = [1, 2, 0.3];
    lb = [0, 0, 0];
    ub = [2, 10000, 1];
    % fitting
    [params_fit,resnorm,residual,exitflag,output,lambda,jacobian] = ...
        lsqnonlin(decay_model, params0, lb, ub);
    

    
    % calculate the decay constant (2nd parameter) and its STD using the
    % 95% confidence interval
    CI = nlparci(params_fit,residual,'jacobian',jacobian);
    decay_constant_STD = (CI(2,2) - CI(2,1))/2;
    data_compiled(i).decay_const = params_fit(2);
    data_compiled(i).decay_const_STD = decay_constant_STD;
    
    hold on
    plot(Time, Fluo)
    plot(Time, params_fit(1)*exp(-params_fit(2)*Time) + params_fit(3))  
%     pause
    
end

%% generate plots for the decay constant (k)
wavelength = extractfield(data_compiled,'wavelength');
decay_const = extractfield(data_compiled,'decay_const');
decay_const_STD = extractfield(data_compiled,'decay_const_STD');
errorbar(wavelength, decay_const, decay_const_STD)

xlabel('wavelength (nm)')
ylabel('export rate (1/min)')
xlim([450 550])
xticks([450 475 500 525 550])
yticks([0 1 2 3 4 5])

box on
StandardFigure(gcf, gca)

saveas(gcf, [FigPath, filesep, 'NLS-mCherry-LEXY-lambda_export_rate.pdf'])
%% Save the result into the structure, data_compiled
save([DataPath 'wavelength_dependence.mat'],'data_compiled')