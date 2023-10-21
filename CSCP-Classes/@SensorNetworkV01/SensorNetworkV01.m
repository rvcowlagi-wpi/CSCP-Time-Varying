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
Class definition of sensor network:
	* Grid locations
	* Placement with basis ID matching
	* Works with a threat model defined by ParametricThreat class

*** AT SOME POINT WE NEED TO FIGURE OUT HOW TO MAKE THIS MORE GENERAL,
E.G., MOST ATTRIBUTES ARE COMMON, ONLY THE SENSOR PLACEMENT/CONFIGURATION
TECHNIQUE WILL DIFFER ***
%}

classdef SensorNetworkV01
	properties
		nSensors
		noiseVariance

		configuration	= [];	% sensor grid locations
		configHistory	= [];	% grid location history
%         allConfigurations = []; % stores all possible configurations
        truepathCost;
        pathLength;
        truepathCostHistory = [];
        estimatedpathCost;
        estimatedpathCostHistory = [];
        varpathCost;
        varpathCostHistory = [];
        pathRisk;
        pathRiskHistory = [];
		identifiedBasis	= [];

		gridWorld		= [];
		threatModel		= [];
        
		% do not use ".threatState" because threatModel is set at
		% initialization and not updated later
	end

	methods
		%==================================================================
		function obj = SensorNetworkV01(nSensors_, noiseVariance_, ...
				threatModel_, gridWorld_)
			% Initialization, including first configuration

            obj.truepathCost = 0;
            obj.truepathCostHistory = [];
            obj.estimatedpathCost = 0;
            obj.estimatedpathCostHistory = [];

			% Set number of sensors
			obj.nSensors		= nSensors_;
			obj.gridWorld		= gridWorld_;
			obj.noiseVariance	= noiseVariance_;
			obj.threatModel		= threatModel_;

			% Random selection of Ns basis from Np^2 basis
			if (obj.nSensors > threatModel_.nStates)
    			randBasis = randperm(threatModel_.nStates, threatModel_.nStates);
			else
    			randBasis = randperm(threatModel_.nStates, obj.nSensors);
			end
			initLocOnGrid = zeros(length(randBasis), 1);

			for m1 = 1:length(randBasis)
    			% Convert basis center to nearest grid point
    			x1_		= threatModel_.basisCenter(1, randBasis(m1));
    			x2_		= threatModel_.basisCenter(2, randBasis(m1));
    			gridVec	= kron( [x1_, x2_]', ones(1, obj.gridWorld.nPoints) ) - ...
					obj.gridWorld.coordinates;

    			x12Sep	= sqrt( gridVec(1,:).^2 + gridVec(2,:).^2 );
    			[~, basisIndex] = min(x12Sep);
    			initLocOnGrid(m1) = basisIndex;
			end

			obj.configuration	= initLocOnGrid;
			obj.configHistory	= initLocOnGrid;

			% Update identified basis
			obj = obj.identify_near_basis();
		end
		%------------------------------------------------------------------
	
%         %==================================================================
%         function observationH_ = calc_rbf_value(obj, configurations_)
% 
% 			% "locations_" is either:
% 			%	2 x n vector, where each column has 2D position coordinates
% 			%	OR
% 			%	n x n x 2 array with meshgrid locations for 2D position	coords
%             locations_	 = obj.gridWorld.coordinates(:, configurations_); 
%             size(locations_);
%             
% 			if size(locations_, 3) > 1
% 				nLocations		= numel(locations_(:, :, 1));
% 				tmpX			= locations_(:, :, 1);
% 				tmpY			= locations_(:, :, 2);
% 				locationsFlatnd = [tmpX(:) tmpY(:)]';
% 			else
% 				nLocations		= size(locations_, 2);
% 				locationsFlatnd	= locations_;
% 			end
% 			observationH_		= zeros(nLocations, obj.threatModel.nStates);
% 
% 		    for m1 = 1:nLocations					
% 				locationVec_ = [locationsFlatnd(1, m1)*ones(1, obj.threatModel.nStates); ...
% 					locationsFlatnd(2, m1)*ones(1, obj.threatModel.nStates)];
% 			
% 				observationH_(m1, :) = exp((-1 / (2 * obj.threatModel.basisSpread)) .* ...
% 					((locationVec_(1, :) - obj.threatModel.basisCenter(1, :)).^2 + ...
% 	 				(locationVec_(2, :) - obj.threatModel.basisCenter(2, :)).^2) );
%                 
%             end
% 		end
% 		%------------------------------------------------------------------

		%==================================================================
		obj = configure(obj,threat_, grid_, optimalPath, timestep)
		% Sensor configuration implemented in a different file
		%------------------------------------------------------------------

        %==================================================================
		obj = configure1(obj,threat_,grid_, optimalPath_, timestep)
		% Sensor configuration implemented in a different file
		%------------------------------------------------------------------

		%==================================================================
		function obj = identify_near_basis(obj)

			% Find identifiable basis functions
			basisIndices = false(1, obj.threatModel.nStates);
			
			% Identifiability defined as within area of significant support: 3*sigma
			maxSep	= 3*sqrt(obj.threatModel.basisSpread);
			for m1 = 1:obj.threatModel.nStates
    			basisVec = kron( ...
					obj.threatModel.basisCenter(:, m1), ...
					ones(1, numel(obj.configHistory)) ) - ...
					obj.gridWorld.coordinates(:, obj.configHistory);

    			% Euclidean distance of basis from sensor locations
    			basisSep = sqrt(basisVec(1,:).^2 + basisVec(2,:).^2);
    			
    			% If the closest sensor to basis is less than 3*sigma, it's identifiable
    			if ( min(basisSep) <= maxSep )
        			basisIndices(m1) = true;
    			end
			end
			tmp00				= 1:obj.threatModel.nStates;
			identifiableBases	= tmp00(basisIndices);
			nBasisIndices		= sum(basisIndices);
			
			%------ Rank basis by closeness to sensors
			% for m1 = 1:numel(identifiable_bases)
			%     basis_sensor_vec = kron( threatModel.basisCenter(:, m1), ...
			% 		ones(1, numel(obj.locationsHistory)) ) - grid_world.coordinates(:, obj.locationsHistory);
			%     basis_sensor_sep(m1) = min(sqrt(basis_sensor_vec(1,:).^2 + basis_sensor_vec(2,:).^2));
			% end
			
			% Find closest basis to each sensor
			for m1 = 1:numel(obj.configHistory)
    			sensorBases = kron( ...
					obj.gridWorld.coordinates(:, obj.configHistory(m1)) , ...
					ones(1, nBasisIndices) ) - ...
					obj.threatModel.basisCenter(:, identifiableBases);

    			[~, nearestBases(m1)] = min(sqrt(sensorBases(1,:).^2 + sensorBases(2,:).^2));
			end
			obj.identifiedBasis = unique(identifiableBases(nearestBases),'stable');
		end
		%------------------------------------------------------------------

        %==================================================================
		obj = plotCost_(obj, flags_)
		% State and estimate plots in a different file
		%------------------------------------------------------------------
	end
end

