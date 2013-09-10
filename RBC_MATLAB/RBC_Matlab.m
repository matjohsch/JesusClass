%% Basic RBC model with full depreciation
%
% Jesus Fernandez-Villaverde
% Haverford, July 31, 2013

%% 0. Housekeeping

clear all
close all
clc

tic

%%  1. Calibration

aalpha = 1/3;     % Elasticity of output w.r.t. capital
bbeta  = 0.95;    % Discount factor
ddelta = 0.09;
param1=aalpha*bbeta/(1-bbeta+ddelta*bbeta);
labourSteadyState =1/3;
pphi = (1-aalpha)/(labourSteadyState^2*(1-ddelta*param1));

% Productivity values
vProductivity = [0.9792; 0.9896; 1.0000; 1.0106; 1.0212]';

% Transition matrix
mTransition   = [0.9727, 0.0273, 0.0000, 0.0000, 0.0000;
                 0.0041, 0.9806, 0.0153, 0.0000, 0.0000;
                 0.0000, 0.0082, 0.9837, 0.0082, 0.0000;
                 0.0000, 0.0000, 0.0153, 0.9806, 0.0041;
                 0.0000, 0.0000, 0.0000, 0.0273, 0.9727];

%% 2. Steady State

labourSteadyState = ((1-aalpha)/(pphi*(1-ddelta*param1)))^(1/2);
capitalSteadyState = (param1)^(1/(1-aalpha))*labourSteadyState;
outputSteadyState = capitalSteadyState^aalpha*labourSteadyState^(1-aalpha);
consumptionSteadyState = outputSteadyState-ddelta*capitalSteadyState;

fprintf(' Output = %2.6f, Capital = %2.6f,Labour = %2.6f, Consumption = %2.6f\n', outputSteadyState, capitalSteadyState, labourSteadyState, consumptionSteadyState); 
fprintf('\n')

% We generate the grid of capital
vGridCapital = 0.8*capitalSteadyState:0.01:1.2*capitalSteadyState;
vGridLabour = 0.8*labourSteadyState:0.01:1.2*labourSteadyState;

nGridCapital = length(vGridCapital);
nGridLabour = length(vGridLabour);
nGridProductivity = length(vProductivity);

%% 3. Required matrices and vectors

%mOutput           = zeros(nGridCapital,nGridLabour,nGridProductivity);
mValueFunction    = zeros(nGridCapital,nGridProductivity);
mValueFunctionNew = zeros(nGridCapital,nGridProductivity);
mPolicyFunctionCapital = zeros(nGridCapital,nGridProductivity);
mPolicyFunctionLabour = zeros(nGridCapital,nGridProductivity);
expectedValueFunction = zeros(nGridCapital,nGridProductivity);

%% 4. We pre-build output for each point in the grid

mCapLab= (vGridCapital'.^aalpha)*(vGridLabour.^(1-aalpha));
mOutput = cat(3,mCapLab*vProductivity(1),mCapLab*vProductivity(2),mCapLab*vProductivity(3),mCapLab*vProductivity(4),mCapLab*vProductivity(5));

%% 5. Main iteration

maxDifference = 10.0;
tolerance = 0.0000001;
iteration = 0;

valueProvisional=zeros(nGridCapital,nGridLabour);
while (maxDifference>tolerance)  
    
    expectedValueFunction = mValueFunction*mTransition';
    
    for nProductivity = 1:nGridProductivity 
        
        for nCapital = 1:nGridCapital
            
            for nCapitalNextPeriod = 1:nGridCapital
                
                for nLabour = 1:nGridLabour
                    
                    consumption = (1-ddelta)*vGridCapital(nCapital)+mOutput(nCapital,nLabour,nProductivity)-vGridCapital(nCapitalNextPeriod);
                    
                    if (consumption<0)
                        valueProvisional(nCapitalNextPeriod,nLabour) = -10000;                  
                    else
                        valueProvisional(nCapitalNextPeriod,nLabour)= log(consumption)-pphi/2*vGridLabour(nLabour)^2+bbeta*expectedValueFunction(nCapitalNextPeriod,nProductivity);
                    end
                end                                  
            end
            [maxV,ind] = max(valueProvisional(:));
            [CapitalChoice,LabourChoice] = ind2sub(size(valueProvisional),ind);
            mValueFunctionNew(nCapital,nProductivity) = maxV;
            mPolicyFunctionCapital(nCapital,nProductivity) = vGridCapital(CapitalChoice);
            mPolicyFunctionLabour(nCapital,nProductivity) = vGridLabour(LabourChoice);
        end
        
    end
    
    maxDifference = max(max(abs(mValueFunctionNew-mValueFunction)));
    mValueFunction = mValueFunctionNew;
    
    iteration = iteration+1;
    if (mod(iteration,10)==0 || iteration ==1)
        fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, maxDifference); 
    end
           
end

fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iteration, maxDifference); 
fprintf('\n')

c=round(length(vGridCapital)/2);
d=round(length(vGridLabour)/2);
fprintf(' My chek = %2.6f\n', mPolicyFunctionCapital(c,3)); 
fprintf('\n')
fprintf(' My chek = %2.6f\n', mPolicyFunctionLabour(d,3)); 
fprintf('\n')
toc

%% 6. Plotting results

figure(1)

subplot(3,1,1)
plot(vGridCapital,mValueFunction)
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
title('Value Function')

subplot(3,1,2)
plot(vGridCapital,mPolicyFunctionCapital)
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
title('Policy FunctionCapital')

subplot(3,1,3)
plot(vGridCapital,mPolicyFunctionLabour)
xlim([vGridCapital(1) vGridCapital(nGridCapital)])
title('Policy FunctionLabour')
