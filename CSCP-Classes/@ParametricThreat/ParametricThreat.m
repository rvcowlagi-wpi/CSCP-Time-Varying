%{
SOFTWARE LICENSE
----------------
Copyright (c) 2023 by 
	Raghvendra V Cowlagi
	Bejamin Cooper
	Prakash Poudel

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
Class definition of parametric threat.
%}

classdef ParametricThreat
	properties
		offset
		alfa
		nStates				% Number of parameters (states)
		basisCenter
		basisSpread

		% For a time-varying threat, the threatState property should be updated
		% at every time step by the "dynamics_discrete" method 
		state

		stateEstimate		% mean state estimate
		estimateCovarPxx	% estimation error covariance
        traceCovarPxx       % trace of error covariance
        pNext                       % covariance prediction for next step
        pNextHist     % covariance prediction for all path_length iterations
        pathLength

		% Maintain histories of state and stateEstimate evolution
		stateHistory
		stateEstimateHistory		% mean state estimate
		estimateCovarPxxHistory		% estimation error covariance
        traceCovarPxxHistory        % trace of error covariance
              
		% Maintain time stamps of state and estimate
		timeStampState
		timeStampEstimate

		% Process noise
		noiseCovarQ

		% Link to a grid world
		ACEGridWorld_

        % State transition matrix
        A
       
	end

	methods
		%==================================================================
		function obj = ParametricThreat(nStates_, halfWorkspaceSize_, ...
				sensorNoiseVar_, grid_)
			% Initialization

			obj.nStates		= nStates_;

			nCenterRows		= sqrt(obj.nStates);
			centerSpacing	= 2 * halfWorkspaceSize_ / (nCenterRows - 1);
			obj.basisCenter	= zeros(2, obj.nStates);
			for m1 = 1:nCenterRows
				for m2 = 1:nCenterRows
					obj.basisCenter(:, (m2 - 1)*nCenterRows + m1) = ...
						[-halfWorkspaceSize_ + (m1 - 1)*centerSpacing; ...
						-halfWorkspaceSize_ + (m2 - 1)*centerSpacing];
				end
			end

			% basis variance chosen such that the range of influence of each basis just
			% touch at the diagonal between two basis centers. Basis centers are 1.4 *
			% center_spacing distance along diagonal, therefore meet at 0.7 * center
			% Range of influce (x_rng) is governed by threat value and noise level. It
			% has similar meaning/use to a signal-to-noise ratio (SNR)
			influenceRange	= 1/2 * sqrt(2) * centerSpacing;
			thetaGuess		= 5;				% theta_guess should be a value in the median range of parameter values
			basisSpread_	= influenceRange^2 / (2 * log(thetaGuess / sensorNoiseVar_));
			obj.basisSpread	= basisSpread_;		% This is \sigma^2_\Psi
			obj.offset		= 1;
% % 
			exampleState	= [ 4    3    4   5   1; ...           % How the threats look visually
								2    5    4   3   5; ...           % but, need to flip, transpose
								6    5    3   5   4;
                              4    2    4   2   3;
                              1    5    3   5   4];
%             exampleState	= [ 4    3    4   5   1   4   1; ...           % How the threats look visually
% 								2    5    4   3   5   2   4; ...           % but, need to flip, transpose
% 								6    5    3   5   4   3   3;
%                               4    2    6   2   3   2   4;
%                               6    5    3   5   4   3   5;
%                               4    2    6   2   3   5   3;
%                               1    5    3   5   4   3   4];
%               exampleState	= [ 4    3    4   5   1    1; ...           % How the threats look visually
% 								2    5    4   3   5    4; ...           % but, need to flip, transpose
% 								6    5    3   5   4    3;
%                               4    2    6   2   3    4;
%                               6    5    3   5   4    5;                          
%                               1    5    3   5   4    4];
%             exampleState	=   [ 4    3    4   5   1   4   1   3; ...           % How the threats look visually
% 								  2    5    4   3   5   2   4   3; ...           % but, need to flip, transpose
% 								  6    5    3   5   4   3   3   2;
%                                   4    2    6   2   3   2   4   1;
%                                   6    5    3   5   4   3   5   1;                                 
%                                   1    5    3   5   4   3   4   3;
%                                   3    2    3   5   4   2   3   1;
%                                   2    2    4   2   3   1   4   4];
%             exampleState	= [ 4    3    4   5   1   4   1  4  3; ...           % How the threats look visually
% 								2    5    4   3   5   2   4  2  3; ...           % but, need to flip, transpose
% 								6    5    3   5   4   3   3  1  2;
%                                 4    2    6   2   3   2   4  6  1;
%                                 6    5    3   5   4   3   5  3  1;
%                                 4    2    6   2   3   5   3  4  2;
%                                 1    5    3   5   4   3   4  1  3;
%                                 3    2    3   5   4   2   3  2  1;
%                                 2    2    4   2   3   1   4  1  4];
%             exampleState	= [ 3    6   1; ...           % How the threats look visually
% 								3    5    2; ...           % but, need to flip, transpose
%                                 1    3    3];
%              exampleState	= [ 3    6   3   1; ...           % How the threats look visually
% 								4    5   2   4; ...           % but, need to flip, transpose
%                                 3    2    4   3
%                                 1    5   2   4];
			obj.state		= reshape(flipud(exampleState)', [obj.nStates, 1]);
			obj.alfa		= 1e-3;


			obj.noiseCovarQ = 0.01*diag( ones(nStates_, 1) );
			obj.ACEGridWorld_= grid_;

			obj.stateEstimate	 = zeros(nStates_, 1);
			obj.estimateCovarPxx = 100 * eye(nStates_);
            obj.traceCovarPxx    = zeros(1,1);

			obj.stateHistory			 = obj.state;
			obj.stateEstimateHistory	 = obj.stateEstimate;
			obj.estimateCovarPxxHistory  = reshape(obj.estimateCovarPxx, nStates_^2, 1);
            obj.traceCovarPxxHistory     = obj.traceCovarPxx;

			obj.timeStampState			= 0;
			obj.timeStampEstimate		= 0;
            obj.pNext                   = 100 * eye(nStates_,nStates_);
            obj.pNextHist               = repmat(100 * eye(nStates_),1,1,obj.ACEGridWorld_.nGridRow^2-1);
%              obj.pathLength              = 9;
%             obj.observationMatrix       = zeros(1,nStates_);
		end
		%---------------------------------------
        % ---------------------------

		%==================================================================
		function c_ = calculate_at_locations(obj, locations_, thisState)
			% "locations_" is either:
			%	2 x n vector, where each column has 2D position coordinates
			%	OR
			%	n x n x 2 array with meshgrid locations for 2D position	coords

			if ~exist("thisState", "var")
				thisState = obj.state;
			end

			observationH= obj.calc_rbf_value(locations_);
			c_			= obj.offset +  observationH * thisState;

			if size(locations_, 3) > 1
				c_		= reshape(c_, size(locations_, 1), size(locations_, 2));
			end

		end
		%------------------------------------------------------------------

		%==================================================================
		function obj = dynamics_discrete(obj, time_step_)

			% This will update the internal state and history. If you just
			% need a prediction for the next time step without changing the
			% state stored in this object, use "process_model" instead.

			t_	= obj.timeStampState(end);

			Phi = obj.calc_rbf_value(obj.basisCenter);
			Ac	= obj.alfa * Phi' ./ (norm(Phi)^2) * ...
				obj.calc_laplacian(obj.basisCenter);
			obj.A	= eye(obj.nStates) + Ac * time_step_;
			
			obj.state = obj.A*obj.state;

			obj.stateHistory	= [obj.stateHistory		obj.state];
			obj.timeStampState	= [obj.timeStampState	t_ + time_step_];
		end
		%------------------------------------------------------------------

		%==================================================================
		function nextState_ = process_model(obj, ...
				threatStatex_, processNoise_, time_step_)
			
			Phi = obj.calc_rbf_value(obj.basisCenter);
			Ac	= obj.alfa * Phi' ./ (norm(Phi)^2) * ...
				obj.calc_laplacian(obj.basisCenter);
			obj.A	= eye(obj.nStates) + Ac * time_step_;						% \theta_dot = A\theta

			nextState_	= obj.A*threatStatex_ + processNoise_;					% Then add noise
		end
		%------------------------------------------------------------------

		%==================================================================
		function c_ = measurement_model(obj, threatStatex_, measNoise_, sensors_)
			% Calculate at locations, then add noise

			locations_	 = obj.ACEGridWorld_.coordinates(:, sensors_.configuration);
			
			observationH = obj.calc_rbf_value(locations_);
			c_			 = obj.offset +  observationH * threatStatex_ + measNoise_;
		end
		%------------------------------------------------------------------
        
		%==================================================================
		function observationH_ = calc_rbf_value(obj, locations_)

			% "locations_" is either:
			%	2 x n vector, where each column has 2D position coordinates
			%	OR
			%	n x n x 2 array with meshgrid locations for 2D position	coords
			if size(locations_, 3) > 1
				nLocations		= numel(locations_(:, :, 1));
				tmpX			= locations_(:, :, 1);
				tmpY			= locations_(:, :, 2);
				locationsFlatnd = [tmpX(:) tmpY(:)]';
			else
				nLocations		= size(locations_, 2);
				locationsFlatnd	= locations_;
			end
			observationH_		= zeros(nLocations, obj.nStates);
%             size(observationH_)

			for m1 = 1:nLocations					
				locationVec_ = [locationsFlatnd(1, m1)*ones(1, obj.nStates); ...
					locationsFlatnd(2, m1)*ones(1, obj.nStates)];
			
				observationH_(m1, :) = exp((-1 / (2 * obj.basisSpread)) .* ...
					((locationVec_(1, :) - obj.basisCenter(1, :)).^2 + ...
	 				(locationVec_(2, :) - obj.basisCenter(2, :)).^2) );
                
			end
        
        end
        
		%------------------------------------------------------------------

		%==================================================================
		function laplacianL_ = calc_laplacian(obj, locations_)
			% Laplacian (d^2/dx^2 + d^2/dy^2)(Gauss basis)

			nLocations	= size(locations_, 2);
			laplacianL_ = zeros(nLocations, obj.nStates);

			for m1 = 1:size(locations_, 2)
				locationVec_ = [locations_(1, m1)*ones(1, obj.nStates); ...
					locations_(2, m1)*ones(1, obj.nStates)];
		
				laplacianL_(m1, :) = (1 / obj.basisSpread) .* ...
					(-2 + 1/(2 * (obj.basisSpread)^2) .* ...
					( (locationVec_(1, :) - obj.basisCenter(1, :)).^2 + ...
					(locationVec_(2, :) -  obj.basisCenter(2, :)).^2))...
        			.* exp((-1 / (2 * (obj.basisSpread)^2)) .* ...
					((locationVec_(1, :) - obj.basisCenter(1, :)).^2 + ...
					(locationVec_(2, :) -  obj.basisCenter(2, :)).^2) );
			end
		end
		%------------------------------------------------------------------

		%==================================================================
		obj = estimate_state_UKF(obj, time_step_, measurementz_k, sensors_)
		% Estimator in a separate file
		%------------------------------------------------------------------

       

        %==================================================================
		obj = estimate_state_UKF1(obj, time_step_, measurementz_k, sensors_, optimalPath)
		% Estimator in a separate file
		%------------------------------------------------------------------

		%==================================================================
		obj = plot_(obj,grid_, flags_)
		% State and estimate plots in a different file
		%------------------------------------------------------------------
	end
end