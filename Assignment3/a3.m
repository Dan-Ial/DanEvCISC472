% Use sphere fitting to calibrate pointer 8700340

stylusID = '8700340';
markerID = '8700339';
dataFile = 'pivot_calibration_0.csv';

setupDrawing();

% Read the raw data into 'pos' (a translation 3-vector) and 'orient'
% (a quaternion 4-vector).

[pos, orient] = read_NDI_data( dataFile, stylusID );

% The following two lines remove the outlier points from the provided
% dataset.  They are used to test the algorithm on "good" data that
% has no outliers.  They should not be used in the final code.

% pos = pos( 200:size(pos,1), : );
% orient = orient( 200:size(orient,1), : );

% fit the sphere to get centre c and radius r

[c, r] = fitSphere( pos );

% fit the sphere using RANSAC, instead.  Note that RANSAC should
% return a list of indices of the inlying poses, and only those poses
% should be used subsequently.

% [c, r, inlierIndices] = fitSphereWithRANSAC( pos );
% pos = pos( inlierIndices, : );
% orient = orient( inlierIndices, : );

% Show the fit

drawCoordSystems( pos, orient );
%drawSphere( c, r );

% Transform c into the coordinate system of each pose

% [YOUR CODE HERE (part 3)]

% Find the average transformed c, which should be the same in all of
% the stylus coordinate systems.  Also find the standard deviation.

% [YOUR CODE HERE (part 3)]

% Report the results
%
% 'c_average' is the average tip position in the stylus coordinate
% system.  'c_stdev' is its standard deviation in the stylus
% coordinate system.

c_average = zeros(1,3);
c_stdev = zeros(1,3);

disp( sprintf( 'Local tip position: (%g, %g, %g)', c_average(1), c_average(2), c_average(3) ) );
disp( sprintf( 'Local tip stdev:    (%g, %g, %g)', c_stdev(1),   c_stdev(2),   c_stdev(3) ) );

% Show vectors to tips in world coordinate system
%
% This is for debugging, so that you can see that the vector touches
% the same pivot point from all stylus coordinate systems.

drawLocalVectorInCoordSystems( c_average, pos, orient );

% Show tip points in global system, along with 95% confidence
% interval as an ellipsoid.
%
% 'c_world' are the tip points in the world coordinate system.
% They should all be very near the pivot point.

c_world = []
drawPointsWithEllipsoid( c_world, c_stdev );



% ---------------- END OF MAIN CODE ----------------



% Fit a sphere to a set of positions
%
% See watkins.cs.queensu.ca/~jstewart/472/notes/15/15-least-squares-problems.html

function [c, r] = fitSphere( pos )

    % [YOUR CODE HERE (part 1)]

    c = zeros(3,1);
    r = 1;
end
  

% Fit a sphere to a set of positions using RANSAC.
%
% ALSO RETURN THE INDICES OF THE BEST INLIERS.  THE CALLING CODE
% SHOULD RESTRICT ITSELF TO THOSE INLIERS.
%
% See en.wikipedia.org/wiki/Random_sample_consensus

function [c, r, bestInlierIndices] = fitSphereWithRANSAC( pos )
    
    % [YOUR CODE HERE (part 2)]

    c = zeros(3,1);
    r = 1;
    bestInlierIndices = 1:size(pos,1);
end
  

% From www.mathworks.com/matlabcentral/fileexchange/35475-quaternions
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


% Read a CSV file from the NDI tracker software and extract the
% position and orientation columns of the named marker.
%
% pos    = n x 3 of translations
% orient = n x 4 of quaternion orientations

function [pos, quat] = read_NDI_data( dataFile, markerID )

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

  pos = zeros( n, 3 );
  quat = zeros( n, 4 );
  k = 1;

  for i = 1:size(t,1)
    if strcmp( t{i,colIndex+3}, 'OK' )

      % Coerce the columns to 'double', since 'readtable' sometimes
      % records numbers as strings.  MATLAB BUG!

      % Extract pose's rotation as a quaternion.  This is in columns
      % offset by +4, +5, +6, +7 from the ID column.

      for j=1:4
	if iscell( t{i,colIndex+3+j} )
	  quat(k,j) = str2double( t{i,colIndex+3+j}{1} );
	else
	  quat(k,j) = t{i,colIndex+3+j};
	end
      end

      % Extract pose's translation as a vector.  This is in columns
      % offset by +8, +9, +10 from the ID column.

      for j=1:3
	if iscell( t{i,colIndex+7+j} )
	  pos(k,j) = str2double( t{i,colIndex+7+j}{1} );
	else
	  pos(k,j) = t{i,colIndex+7+j};
	end
      end

      k = k + 1;
    end
  end
  
  disp( sprintf( '%d points collected', size(pos,1) ) );
end



% Set up the drawing

function setupDrawing() 

    f = figure(1);
    clf(f);
    view(3);
    daspect( [1 1 1] );
    pbaspect manual;
    hold on;
end


% Draw a set of coordinate systems
%
% pos and orient store their vectors in rows.

function drawCoordSystems( pos, orient )
    
    colours = [ 'r' 'g' 'b' ];

    scale = 0.005 * norm(max(pos) - min(pos));

    for i=1:size(pos,1)
        m = quaternion_to_matrix( orient(i,:) );
        t = pos(i,:);
        for j=1:3
            head = t + m(:,j)' .* scale;
            plot3( [t(1) head(1)], [t(2) head(2)], [t(3) head(3)], colours(j) );
        end
    end
end


% Draw a sphere

function drawSphere( c, r )
    
    [x,y,z] = sphere;
    surf( x*r+c(1), y*r+c(2), z*r+c(3), 'FaceAlpha', 0.05, 'FaceColor', [0.6 0.3 0.3] );
end


% Draw a local vector in different coordinate systems
%
% v is a row vector.  pos and orient store their vectors in rows, too.

function drawLocalVectorInCoordSystems( v, pos, orient )
    
    for i=1:size(pos,1)
        m = quaternion_to_matrix( orient(i,:) );
        t = pos(i,:);
        for j=1:3
            head = t + (m * v')';
            plot3( [t(1) head(1)], [t(2) head(2)], [t(3) head(3)] );
        end
    end
end


% Draw points with a 95% CI ellipsoid.
%
% Use matlab's 'ellipsoid' and 'surf' functions.

function drawPointsWithEllipsoid( points, stdev )
    
    % [YOUR CODE HERE (part 4)]

end
