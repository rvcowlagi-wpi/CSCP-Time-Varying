%{
SOFTWARE LICENSE
----------------
Copyright (c) 2023 by Raghvendra V. Cowlagi

Permission is hereby granted to any person obtaining a copy of this
software and associated documentation files (the "Software"), to deal in
the Software, including the rights to use, copy, modify, merge, copies of
the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:  

* The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.
* The Software, and its copies or modifications, may not be distributed,
published, or sold for profit. 
* The Software, and any substantial portion thereof, may not be copied or
modified for commercial or for-profit use.

The software is provided "as is", without warranty of any kind, express or
implied, including but not limited to the warranties of merchantability,
fitness for a particular purpose and noninfringement. In no event shall the
authors or copyright holders be liable for any claim, damages or other
liability, whether in an action of contract, tort or otherwise, arising
from, out of or in connection with the software or the use or other
dealings in the software.      


PROGRAM DESCRIPTION
-------------------
A modular implementation of coupled sensor configuration and planning
(CSCP).
%}


%% Tabula Rasa
clear variables; close all; clc
addpath(genpath('CSCP-Classes'))

%% Initialize

%----- Time and iteration counters
% Time counter: there can be multiple and possibly variable time steps per
% iteration
time_k	= 0;	
k		= 0;	% Iteration counter

%----- Problem dimensions
N_THREAT_STATE	= 9;			
N_SENSORS		= 1;
N_GRID_ROW		= 5;

%----- Other
N_EXP_ITER		= 1;	% Ballpark of how many CSCP iterations may be needed
N_MAX_ITER		= 10;	% Terminate if this is exceeded
TERM_PLAN_RISK	= 0.1;	% Risk threshold for terminating CSCP iterations
SENSOR_NOISE_VAR= 0.1;	% Variance of (i.i.d.) measurement noise in each sensor, assuming homogeneous sensors

time_step_		= 0.1; % ** THIS MAY CHANGE DURING THE LOOP; FIX LATER

%----- Instantiate sensor and threat classes
grid_			= ACEGridWorld(1, N_GRID_ROW);
threat_			= ParametricThreat(N_THREAT_STATE, ...
	grid_.halfWorkspaceSize, SENSOR_NOISE_VAR, grid_);
% threat_         = threat_.estimate_state_UKF1(time_step_,);
sensor_		    = SensorNetworkV01(N_SENSORS, ...
	             SENSOR_NOISE_VAR, threat_, grid_);

%----- Storage for results at each iteration
timeStampMeas	= zeros(1, N_EXP_ITER);
measurementz	= zeros(N_SENSORS, N_EXP_ITER);
planState		= zeros(N_GRID_ROW^2, N_EXP_ITER);		% Planned path
planCostRisk	= zeros(2, N_EXP_ITER);					% Expected cost and risk



%% CSCP Loop
while (1)
	%----- Increment iteration counter
	k	   = k + 1;
	time_k = time_k + time_step_;

	%----- Configure sensors
	sensor_			= sensor_.configure([]);
   
    %----- Propagate true threat
	threat_			= threat_.dynamics_discrete(time_step_);
	trueThreat_k	= threat_.calculate_at_locations( ...
		grid_.coordinates(:, sensor_.configuration) );

    %----- Simulate sensor measurements
	measNoise_k		= sqrt(SENSOR_NOISE_VAR) * randn( N_SENSORS, 1);
	measurementz_k	= trueThreat_k + measNoise_k;	% Pointwise measurement of threat

	%----- Run estimator
	threat_			 = threat_.estimate_state_UKF(time_step_, measurementz_k, sensor_);
    threatStateHat_k = threat_.stateEstimate;
	
    %----- Find optimal plan
% 	grid_			= grid_.plan_path(threat_);
% 	planState_k		= grid_.optimalPath;
% 	planCostRisk_k	= [grid_.pathCost; grid_.pathRisk];

	%----- Store results of this iteration
	timeStampMeas(:, k)		= time_k;
	measurementz(:, k)		= measurementz_k;
% 	planState(:, k)			= planState_k;
% 	planCostRisk(:, k)		= planCostRisk_k;
   
    %----- Run estimator
% 	threat_			 = threat_.estimate_state_UKF1(time_step_, sensor_);
	
	%----- Check termination criteration and break
	if (grid_.pathRisk <= TERM_PLAN_RISK) || (k >= N_MAX_ITER), break; end
end

%% Plot results

% flags_.SHOW_TRUE	 = true;
% flags_.SHOW_ESTIMATE = true;
% threatStatePlotAxes  = threat_.plot_(flags_);


flags_.SHOW_TRUE	= true;
flags_.SHOW_ESTIMATE= false;
flags_.JUXTAPOSE	= true;
grid_.plot_parametric(threat_, sensor_, flags_)