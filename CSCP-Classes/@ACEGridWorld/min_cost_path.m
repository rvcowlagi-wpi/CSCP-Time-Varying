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
Path optimization in the class definition of grid world. "Discovers"
neighbors and assigns vertex numbers as they are discovered.
%}

function obj = min_cost_path(obj)


%------ Search settings
isAllGoalSearch	= obj.searchSetup.reachAllGoals;
idStart			= obj.searchSetup.start;
idGoal			= obj.searchSetup.goal;

%------ Cost, heuristic, and neighbor functions
fcn_heuristic	= obj.searchHeuristicFcn;
fcn_nhbr		= obj.neighborFcn;

nVertices		= 1E3;														% For large graphs, initialize search data structure with arbitrary size
knownVertices	= zeros(nVertices, 1);

% The assumption here is that the function that discovers neighboring
% vertices also provides a **unique** identifier to each vertex. The main
% difference from a typical Dijkstra implementation is that the vertex ID
% (which can be large) need not match this algorithm's internal indexing,
% which is helpful when the #vertices searched << total number of verices.
tmp				= struct('id', 0, 'mk', 0, 'd', Inf, 'b', [], 'x', []);
vertexData		= repmat(tmp, 1, nVertices);

% ID 'id' and state 'x' are redundant; 'ID' is scalar

vertexData(1).id= idStart;
vertexData(1).mk= 1;
vertexData(1).d	= 0;

nKnownVertices	= 1;
knownVertices(1)= idStart;

fringe_			= [1 0];
nFringe			= 1;

isGoalClosed	= 0;

nIter	= 0;
while (nFringe ~= 0) && (~isGoalClosed)
% 	clc;
% 	fprintf('Number of iterations  : %i\n', nIter);
% 	fprintf('Number of OPEN vertexs  : %i\n', nOpen);
% 	fprintf('Number of known vertexs : %i\n\n', nKnownvertexs);
	
	nIter		= nIter + 1;
	vCurrent	= fringe_(1, 1);											% Get vertex from top of (sorted) open stack
	vertexData(vCurrent).mk = 2;											% Mark that vertex as dead
	
	nFringe		= nFringe - 1;
	fringe_(1, :) = [];	
														
	[nhbrID, nhbrCosts] = fcn_nhbr(vertexData(vCurrent).id);				% Other function handle for nhbrs and costs
	
	for k = 1:numel(nhbrID)													% For all neighbors

		if any(nhbrID(k) == knownVertices(1:nKnownVertices))
			tmp2		= 1:nKnownVertices;
			vNeighbor	= tmp2(nhbrID(k) == knownVertices(1:nKnownVertices));
			vNew		= vNeighbor;
		else
			vNew		= nKnownVertices + 1;
			vertexData(vNew).id	= nhbrID(k);
			vertexData(vNew).mk	= 0;
			
			nKnownVertices		= vNew;
			knownVertices(vNew)	= nhbrID(k);	
		end
		costNew	= nhbrCosts(k);												% Cost to go from act to new
		
% 		[v_current v_new cost_new]
		
		if vertexData(vNew).mk == 0											% Unvisited
			vertexData(vNew).mk	= 1;										% Mark open
			vertexData(vNew).d	= vertexData(vCurrent).d + costNew;			% Update c2come of newly visited state
			vertexData(vNew).b	= vCurrent;
			

			tmp_open = bin_sort(fringe_(1:nFringe, :), [vNew vertexData(vNew).d], 2);
% 			if numel(tmp_open) == 0
% 				nFringe	= 0;
% 				fringe_	= [];
% 			else
				nFringe	= size(tmp_open, 1);
				fringe_(1:nFringe, :)	= tmp_open;							% Add [vertexNew cost] to sorted open list
% 			end			
		elseif vertexData(vNew).mk == 1										% Already open, update c2come if necessary
			if vertexData(vNew).d > vertexData(vCurrent).d + costNew
				vertexData(vNew).d	= vertexData(vCurrent).d + costNew;
				vertexData(vNew).b	= vCurrent;
				
				fringe_( (fringe_(1:nFringe, 1) == vNew), :) = [];
				nFringe			= nFringe - 1;
				
				tmp_open		= bin_sort(fringe_(1:nFringe, :), ...
					[vNew vertexData(vNew).d], 2);					

% 				if numel(tmp_open) == 0
% 					nFringe		= 0;
% 					fringe_	= [];
% 				else
					nFringe		= size(tmp_open, 1);
					fringe_(1:nFringe, :)	= tmp_open;						% Add [vertexNew cost] to sorted open list
% 				end
			end
		end
	end
	

	if ~isAllGoalSearch
		for k = 1:numel(idGoal)
			[is_goal_known, idx_goal] = ismember(idGoal(k), knownVertices);
			if is_goal_known && vertexData(idx_goal).mk == 2, isGoalClosed = 1; break; end
		end
	else			
		if ~all(ismember(idGoal, knownVertices)), continue; end
		isGoalClosed = 1;
		for k = 1:numel(idGoal)
			[~, idx_goal] = ismember(idGoal(k), knownVertices);
			if vertexData(idx_goal).mk == 2, isGoalClosed = 0; break; end
		end
	end

end

knownVertices = knownVertices(1:nKnownVertices);

%----- Outputs
obj.optimalPath = zeros(obj.nPoints, 1);
obj.optimalPath(1:5) = 1:5;
obj.pathCost	= 1;
obj.pathRisk	= 2;

end

