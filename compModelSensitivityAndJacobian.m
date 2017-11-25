%% info: http://sbml.org/Software/SBMLToolbox/SBMLToolbox_4.0_API_Manual

model='test_model/BIOMD0000000404.xml';
model_id = 'BIOMD0000000404';

addpath(pwd);

if ~exist(model_id, 'dir')
    mkdir(model_id);
end

cd(model_id);

%% Load model
[SBMLModel, errors] = TranslateSBML(['../' model],1,0);


% get the name/id of the model

Name = '';
if (SBMLModel.SBML_level == 1)
    Name = SBMLModel.name;
else
    if (isempty(SBMLModel.id))
        Name = SBMLModel.name;
    else
        Name = SBMLModel.id;
    end;
end;

if (length(Name) > 63)
    Name = Name(1:60);
end;

model_name = Name;
model_handle = [ '@' model_name ];




%% Get variable names
SpeciesNames = GetSpecies(SBMLModel);
[VarParams, VarInitValues] = GetVaryingParameters(SBMLModel);

VarNames = [ SpeciesNames VarParams ];

%% Write ODE function
% Modified version to add a perturbation
% WriteODEFunctionForSensitivity(SBMLModel);
% WriteODEFunction(SBMLModel);

%% Steady-state and Jacobian
% initial concentrations
x0 = eval(model_name);

% go to steady state without perturbation and save it
[t,x] = ode23tb(eval(model_handle), [0, 10], x0);
xss = x(end,:);% + randn() * 1e-15;

% Save steady state
csvwrite('steady_state.csv', xss');

% Compute Jacobian
J = compJacobian(@(x) eval([ model_name '(10, x)']),xss);
Jtab = array2table(J, 'VariableNames', VarNames, 'RowNames', VarNames);
writetable(Jtab, 'jacobian.csv', 'WriteRowNames', true);


cd('..')