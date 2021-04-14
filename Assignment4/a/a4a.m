% Perform ICP to align collected points with a model
%
%
% The 'collectedFilename' contains a stream of points collected from
% the stylus.  But we should use only those points at which the stylus
% stays steady for a while.
%
% The 'femur.stl' is the 3D computer model of the femur that we are
% registering to.


modelFilename = 'femur.stl';            % model
collectedFilename = 'knee2.csv';        % collected points
stylusID = '8700340';                   % ID in collected points file
stylusTip = [-17.02 -1.23 -157.13];     % Calibrated stylus tip position from Assignment 3


% number of attempts running ICP with different initial pose

numAttempts = 1;   % CHANGE THIS AFTER GETTING apply_ICP() WORKING


% Get the stylus-collected points

pts = read_NDI_data( collectedFilename, stylusID, stylusTip );
disp( sprintf( 'read %d stylus points', length(pts)) );

observedPts = find_observed_points( pts ); % only keep those at which the stylus pauses
n = size(observedPts,1);

disp( sprintf( 'found %d points at which the stylus paused', length(observedPts)) );


% read the model and get its unique vertices

[modelFaces, modelPts] = stlread( modelFilename );
modelPts = unique( modelPts, 'rows' ); % remove duplicate points in model and discard faces

disp( sprintf( 'read %d model points', length(modelPts)) );


% Build KD-tree from model vertices

kdTree = KDTreeSearcher( modelPts, 'BucketSize', 10 );


% Run ICP several times on the observedPts

for i = 1:numAttempts

  % Find an initial rotation and translation of the observedPts that
  % are PLAUSIBLE.
  %
  % For example, you could find the principal axes of the modelPts,
  % then translate the observedPts to lie within the 95% confidence
  % interval of the point cloud's mean.
  %
  % Any rotation works, so you could use a *uniform* *random*
  % rotation matrix.
  %
  % Alternatively, you could try rotations that align the principal 
  % axes of the observed points to be NEAR the principal axes of
  % the model points.

  % Pick an initial translation
  % 
  % [YOUR CODE HERE (after you get apply_ICP working)]

  initTrans = [0 0 0];

  % Pick a uniform random rotation
  %
  % [YOUR CODE HERE (after you get apply_ICP working)]

  initRot = eye(3,3);
  
  % Apply ICP with this initRot and initTrans

  disp( sprintf( "\nAttempt %d", i ) );

  [accumRot,accumTrans,rmsError] = apply_ICP( observedPts, initRot, initTrans, kdTree, modelPts )

  % Keep the best so far

  if rmsError < bestRMSError || i == 1

    bestRMSError = rmsError;
    bestRot = accumRot;
    bestTrans = accumTrans;
  end
end


% Report the final simplified transform.  The final transform is
%
%    R * (pts - mu) + mu + t
%
% which simplifies to
%
%    R * pts + (- R * mu + mu + t)

mu = mean(observedPts);

R = bestRot
t = - (R * mu')' + mu + bestTrans
bestRMSError


% Show the best solution

xPts = (R * observedPts')' + t;

indices = knnsearch( kdTree, xPts );
closestPts = modelPts(indices,:);

draw_all( modelPts, xPts, closestPts );


% ---- Done ----


% Run ICP on some points to which an initial rotation and translation
% are applied.  Return the total rotation and translation after ICP
% terminates.  We'll presume that rotations are always about the
% centroid of the points.


function [R,t,rmsError] = apply_ICP( pts, initRot, initTrans, kdTree, modelPts )

  % parameters for one run of ICP

  maxRMSError = 0;             % ICP error (in mm)
  maxIterations = 200;           % ICP iterations
  noImprovementFraction = 0.001; % quit iterating when RMSE improves by less than this fraction
  
  % Apply initRot and initTrans to points.  Note that rotation is
  % applied around the centroid of the points.
  
  xPts = (initRot * (pts - mean(pts))')' + mean(pts) + initTrans;

  % Accumulate a transformation over multiple iterations of ICP in this function.
  %
  % This transformation will transform OBSERVED points onto MODEL points.

  accumRot   = initRot;     % accumulated rotation 
  accumTrans = initTrans;   % accumulated translation
  prevRMSE = 0;

  rmsError = 9999;    % current error
  iter = 0;           % current iteration

  while rmsError > maxRMSError & iter < maxIterations & abs(rmsError - prevRMSE) > noImprovementFraction*prevRMSE
    % Find the nearest model points using kD-tree

    indices = knnsearch( kdTree, xPts );
    closestPts = modelPts(indices,:);

    % [ YOUR CODE HERE ]
    % With each iteration in the ICP loop, you should accumulate the
    % incremental translations and rotations into 'accumTrans' and 'accumRot'
    
    prevRMSE = rmsError;%moved
    
    % 1. for each point (p) in the point set (P), find the closest model 
    % point (m) in the model (M)
    % already (mostly) done above in the given skeleton code
    % P = xPts
    % p = indices
    % M = modelPts
    % m = closestPts
    
    % 2. given the paired points (p,m), apply Procrustes to find  a
    % transformation (T') , to move the p to the corresponding m
    
    [m, n] = size(xPts);
    % inital rotation/translation
    % accumRot and accumTrans
    
    % let Q = closestPts'
    Q = closestPts';
    
    % use Procrustes method to get rotation
    meanP = repmat( mean(xPts',2), 1, m );
    meanQ = repmat( mean(Q,2), 1, m );
    
    Pcentred = xPts' - meanP;
    Qcentred = Q - meanQ;

    [U, S, V] = svd( Qcentred * Pcentred' );

    accumRot2 = U * V';
    
    % caluclate translation
    accumTrans2 = meanQ - accumRot2 * meanP;
    
    % 3. add T' to the accumulated transformation:  T <- T'T
    accumRot = accumRot2' + accumRot;
    accumTrans = accumTrans2' + accumTrans;

    % 4. update points using the transformation: p <- T'p
    xPts = ((accumRot2 * xPts') + accumTrans2)';
    
    % 5. RMS error: determine the error between each p and the corresponding m
    rmsError = sqrt(mean((xPts(:)-closestPts(:)).^2));
    
    % Update the display

    draw_all( modelPts, xPts, closestPts );

    disp( sprintf( '%2d: RMSE = %.2f', iter, rmsError ) );

    iter = iter + 1;
  end % end while

  
  R = accumRot;
  t = accumTrans;

  % Report

  if abs(rmsError - prevRMSE) <= noImprovementFraction*prevRMSE
    disp( sprintf( 'Stopped due to insufficient (%.2f%%) improvement.', abs(rmsError-prevRMSE)/prevRMSE*100 ) );
  elseif iter >= maxIterations
    disp( sprintf( 'Stopped after limit of %d iterations was reached.', maxIterations ) );
  else	
    disp( 'Stopped with low RMS error.' );
  end

end



% In a stream of 3D points, find those that are in approximately the same
% position for a "long time".

function ptsOut = find_observed_points( ptsIn )

  minRestingPts = 20;   % need at least this many observed resting points
  maxRestingDist = 2.0; % max distance (mm) between observed resting points
  
  ptsOut = [];
  pt0 = [0 0 0];
  count = 0;
  
  for k = 1:length(ptsIn) % Note: this will discard the last set of resting points
                          % when it's not followed by a non-resting point
    pt = ptsIn(k,:);
    dist = norm( pt-pt0 );
    if dist <= maxRestingDist % add a point
      count = count+1;  
    elseif count < minRestingPts % end of resting points ... too few
      pt0 = pt;
      count = 0;
    else % end of resting points ... enough
      ptsOut = [ ptsOut; median( ptsIn( k-count:k-1, : ) ) ];
      count = 0;
    end
  end
end


% Draw the model points, the current (transformed) observed points,
% and lines from the observed points to their closest model points.

function draw_all( modelPts, observedPts, closestPts )
    
  f = figure(1);
  clf(f);
  
  axis('vis3d');
    
  % Draw the model points
  
  scatter3( modelPts(:,1), modelPts(:,2), modelPts(:,3), 1, [0.8, 0.8, 0.8] ); % grey

  % Draw the transformed observed points
  
  hold on;
  
  scatter3( observedPts(:,1), observedPts(:,2), observedPts(:,3), 30, [0.8, 0.3, 0.1] ); % reddish
  
  % Draw lines between correspondings points
  
  for i = 1:length(observedPts)
    plot3( [observedPts(i,1), closestPts(i,1)], [observedPts(i,2), closestPts(i,2)], [observedPts(i,3), closestPts(i,3)] ); 
  end
  
  drawnow;
  hold off;
end


% Read a CSV file from the NDI tracker software and extract the
% position and orientation columns of the named marker.
%
% pos    = n x 3 of translations
% orient = n x 4 of quaternion orientations

function pts = read_NDI_data( dataFile, markerID, localVector )

  t = readtable( dataFile, 'PreserveVariableNames', true );

  % Find the column 

  colIndex = find(contains( t.Properties.VariableNames, markerID ));

  if colIndex == []
    disp( sprintf( "In %s: Could not find a column header containing '%s'.", dataFile, markerID ) );
    exit;
  end

  % From the rows with state == 'OK' (state is at +3 offset from the ID
  % column), extract the Tx, Ty, Tz columns.

  status = t{:,colIndex+3};
  n = size( t( strcmp(status,'OK'), 1 ), 1 );

  pts = zeros( n, 3 );
  k = 1;

  for i = 1:size(t,1)
    if strcmp( t{i,colIndex+3}, 'OK' )

      % Coerce the columns to 'double', since 'readtable' sometimes
      % records numbers as strings.  MATLAB BUG!

      % Extract pose's rotation as a quaternion.  This is in columns
      % offset by +4, +5, +6, +7 from the ID column.

      quat = zeros( 1, 4 );
      for j=1:4
        if iscell( t{i,colIndex+3+j} )
          quat(j) = str2double( t{i,colIndex+3+j}{1} );
        else
          quat(j) = t{i,colIndex+3+j};
        end
      end

      % Extract pose's translation as a vector.  This is in columns
      % offset by +8, +9, +10 from the ID column.

      pos = zeros( 1, 3 );
      for j=1:3
        if iscell( t{i,colIndex+7+j} )
          pos(j) = str2double( t{i,colIndex+7+j}{1} );
        else
          pos(j) = t{i,colIndex+7+j};
        end
      end
      
      pts(k,:) = (quaternion_to_matrix( quat ) * localVector')' + pos;

      k = k + 1;
    end
  end
end


% From https://www.mathworks.com/matlabcentral/fileexchange/35475-quaternions
%
% Convert quaternion to 3x3 matrix.

function R = quaternion_to_matrix( Qrotation )

  w = Qrotation( 1 );
  x = Qrotation( 2 );
  y = Qrotation( 3 );
  z = Qrotation( 4 );

  Rxx = 1 - 2*(y^2 + z^2);
  Rxy = 2*(x*y - z*w);
  Rxz = 2*(x*z + y*w);
  Ryx = 2*(x*y + z*w);
  Ryy = 1 - 2*(x^2 + z^2);
  Ryz = 2*(y*z - x*w );
  Rzx = 2*(x*z - y*w );
  Rzy = 2*(y*z + x*w );
  Rzz = 1 - 2 *(x^2 + y^2);

  R = [ Rxx,    Rxy,    Rxz;
        Ryx,    Ryy,    Ryz;
        Rzx,    Rzy,    Rzz  ];
end
