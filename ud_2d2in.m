function [DEFL,REACT,ELE_FOR,ELE_YLD,AFLAG,APRATIOS,LIMIT_STATE] = ud_2d2in(...
				nnodes,coord,concen,fixity,nele,ends,A,Izz,Zzz,Ayy,...
				E,v,Fy,YldSurf,Wt,webdir,w,thermal,truss,E_t,anatype,sol_scheme,numsteps,...
				ratio_req,stop_ratio,restart,defl,react,ele_for,ele_yld,...
				apratios,limit_state,h_stat_mes)
%UD_2D2IN performs a user defined two-dimensional
% second-order inelastic analysis of a structural system.
% It is assumed that all information is defined in the
% X-Y plane.
%

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Functions Called
%              < to be defined by the student >
%
%  Dictionary of Variables
%     Input Information:
%       nnodes         ==  total number of nodes
%       coord(i,1:2)   ==  node i's coordinates
%                            coord(i,1) = X coordinate
%                            coord(i,2) = Y coordinate
%       concen(i,1:3)  ==  concentrated loads for node i's 3 d.o.f.
%                            concen(i,1) = force in global X direction
%                            concen(i,2) = force in global Y direction
%                            concen(i,3) = moment about global Z axis
%       fixity(i,1:3)  ==  prescribed displacements for node i's 3 d.o.f.
%                          Note: A free d.o.f. will have a value of NaN
%                          and hence, you will find the Matlab function
%                          isnan very useful.
%                          Examples: If fixity(15,1) is set to NaN, then node 15's
%                                      X-disp component is free;
%                                    If fixity(2,3) is set to 0.0, then node 2's
%                                      Z-rotation component is supported;
%                                    If fixity(5,2) is set to -2.1, then node 5's
%                                      Y-disp component is supported and defined
%                                      with a settlement of -2.1 units.
%                            fixity(i,1) = prescribed disp. in global X direction
%                            fixity(i,2) = prescribed disp. in global Y direction
%                            fixity(i,3) = prescribed rotation about global Z axis
%       nele           ==  total number of elements
%       ends(i,1:8)   ==  element i's nodal information
%                            ends(i,1) = start node #
%                            ends(i,2) = finish node #
%                            ends(i,3) = flag to indicate whether or not flexural
%                            moment is released at start node.  ends(i,3)=0 not
%                            released (rigid connection); ends(i,3)=1 flexural
%                            moment is released (pinned connection); ends(i,3)=3
%                            at least one of the flexural moments are partially or fully
%                            released (see below for connection stiffness attributes)
%                            ends(i,4) = flag to indicate whether or not flexural
%                            moment is released at finish node.  ends(i,4)=0 not
%                            released (rigid connection); ends(i,4)=1 flexural
%                            moment is released (pinned connection); ends(i,4)=2
%                            at least one of the flexural moments are partially or fully
%                            released (see below for connection stiffness attributes)
%                            ends(i,5) = rotational spring stiffness at the start
%                            node and about element i's local z-z axis.
%                            ends(i,6) = rotational spring stiffness at the finish
%                            node and about element i's local z-z axis.
%                            ends(i,7) = connection moment capacity Mpz at the start
%                            node and about element i's local z-z axis.
%                            ends(i,8) = connection moment capacity Mpz at the finish
%                            node and about element i's local z-z axis.
%       A(i)           ==  element i's cross sectional area
%       Izz(i)         ==  element i's moment of inertia about its local z-z axis
%       Zzz(i)         ==  element i's plastic section modulus about its local z-z axis
%       Ayy(i)         ==  element i's effective shear area along its local y-y axis
%       E(i)           ==  element i's material elastic modulus, Young's Modulus
%       v(i)           ==  element i's material Poisson's ratio
%       Fy(i)          ==  element i's material yield strength
%       YldSurf(i)     ==  element i's yield surface maximum values
%                              YldSurf(i,1) = maximum P/Py value
%                              YldSurf(i,2) = maximum Mz/Mpz value
%                              YldSurf(i,3) = maximum My/Mpy value
%       Wt(i)          ==  element i's material weight density
%       webdir(i,1:3)  ==  element i's unit web vector.  This is a unit vector
%                          that defines the element's local y-y axis with respect
%                          to the global coordinate system.  It is based only on the
%                          structures undeformed geometry.
%                            webdir(i,1) = x component of element's unit web vector
%                            webdir(i,2) = y component of element's unit web vector
%                            webdir(i,3) = z component of element's unit web vector
%                                          (which will always be zero for 2D analysis)
%                          NOTE: An element's 3x3 rotation matrix, [g], is constructed
%                          as follows: First, calculate a unit vector, x_vect, that
%                          describes the element's local x-axis (be sure to include
%                          all three components with the z component always being zero).
%                          Second, take the cross product of x_vect and webdir(i,:) to
%                          obtain z_vect, i.e. z_vect = cross(x_vect,webdir(i,:)). Third,
%                          set z_vect to a unit vector, i.e. z_vect = z_vect/norm(z_vect).
%                          For 2D analysis, z_vect will be either [0 0 1] or [0 0 -1].
%                          Finally, the first row of [g] is x_vect, its second row is
%                          webdir(i,:), and its third row is z_vect.
%       w(i,1:2)       ==  element i's uniform load which references its
%                            local coordinate system
%                              w(i,1) = x component of uniform load
%                              w(i,2) = y component of uniform load
%       thermal(i,1:4)   ==  element i's thermal strain effects which reference its
%                            local coordinate system
%                              thermal(i,1) = coefficient of thermal expansion
%                              thermal(i,2) = change in temperature at centroid
%                              thermal(i,3) = linear temperature gradient in local y-dir
%                                           = (T_up_y - T_btm_y) / depth_y
%                              thermal(i,4) = linear temperature gradient in local z-dir
%                                           = (T_up_z - T_btm_z) / width_z
%       truss          ==  flag to indicate if structure is a truss or not
%                              truss = 0   System is not a truss
%                              truss = 1   System is a truss
%       E_t              ==  flag to indicate if tangent modulus approach should
%                            be employed
%                              E_t = 1   Do NOT employ tangent modulus theory
%                              E_t = 2   Employ tangent modulus theory
%       anatype        ==  flag to indicate which type of analysis is requested
%                              anatype = 1  First-Order Elastic
%                              anatype = 2  Second-Order Elastic
%                              anatype = 3  First-Order Inelastic
%                              anatype = 4  Second-Order Inelastic
%                              anatype = 5  Elastic Buckling (Eigenvalue)
%                              anatype = 6  Inelastic Buckling (Eigenvalue)
%       sol_scheme     ==  flag to indicate which type ofsolution scheme is requested
%                              sol_scheme = 1   Simple Step (Euler)
%                              sol_scheme = 2   Predictor-Corrector (Midpoint R-K)
%                              sol_scheme = 3   Other 1
%                              sol_scheme = 4   Other 2
%                              sol_scheme = 5   Other 3
%       numsteps       ==  requested maximum number of load steps or increments
%       ratio_req      ==  requested load step or increment size
%       stop_ratio     ==  requested maximum applied load ratio
%       restart        ==  flag to indicate if a new or continuing analysis
%                              restart = 1  Start a new analysis
%                              restart = 2  Continue with previous analysis results
%       defl(i,1:3,n)    ==  node i's calculated 3 d.o.f. deflections for all
%                            n previous steps.  See DEFL for additional details.
%       react(i,1:3,n)   ==  node i's calculated 3 d.o.f. reactions for all
%                            n previous steps. See REACT for additional details.
%       ele_for(i,1:6,n) ==  element i's calculated 6 element end forces for all
%                            n previous steps.  See ELE_FOR for additional details.
%       ele_yld(i,1:2,n) ==  applied load ratio at which element i's start and/or end node
%                            formed a plastic hinge for all n previous steps.  See ELE_YLD
%                            for additional details.
%       apratios(i,1)    ==  applied load ratio at previous step i
%                            Note that apratios is a COLUMN vector.
%                            size(apratios,1) = total number of previous load steps.
%                            See APRATIOS for additional details.
%       limit_state      ==  flag indicating post-limit state analysis
%                              limit_state = 0  pre-limit state, system is loading
%                              limit_state = 1  post-limit state, system is unloading
%       h_stat_mes       ==  handle to graphics object which displays status message.
%                            Example:  % A message at the end of each load step is
%                                      % often informative to the user.  The following
%                                      % Matlab code displays a message which
%                                      % includes the contents of variable istep
%                                      text_mess = ['Performing step #',num2str(i_step)];
%                                      set(h_stat_mes,'String',text_mess); drawnow;
%
%     Local Information:
%              < to be defined by the student >
%
%     Output Information:
%       DEFL(i,1:3,m)    ==  node i's calculated 3 d.o.f. deflections for all m steps
%                            of the analysis.  If this is a new analysis, DEFL should be
%                            set equal to [].  If this is a continuing analysis, DEFL
%                            should be first set equal to defl, i.e. DEFL = defl;
%                            For a supported d.o.f., the proper percent of the prescribed
%                            value, fixity(i,dof), should be inserted in DEFL.
%                            Example: DEFL(12,2,5) is node 12's Y-disp. component
%                                     at the end of the fifth load increment.
%                                     If the Y-disp. d.o.f. for node 12
%                                     is supported, then the appropriate percent of
%                                     the prescribed value, fixity(12,2), should
%                                     inserted at DEFL(12,2,5).
%                              DEFL(i,1,m) = displacement in X direction at end of step m
%                              DEFL(i,2,m) = displacement in Y direction at end of step m
%                              DEFL(i,3,m) = rotation about Z direction at end of step m
%       REACT(i,1:3,m)   ==  node i's calculated 3 d.o.f. reactions for all m steps
%                            of the analysis.  If this is a new analysis, REACT should be
%                            set equal to [].  If this is a continuing analysis, REACT
%                            should be first set equal to react, i.e. REACT = react;
%                            For free d.o.f., a value of 0.0 should be inserted in REACT. 
%                            Example: REACT(5,3,8) is node 5's Z-moment reaction
%                                     component at the end of the eighth load
%                                     increment.  If the Z-rotation d.o.f. for
%                                     node 5 is free, then a value of 0.0 should
%                                     inserted at REACT(5,3,8).
%                              REACT(i,1,m) = force in X direction at end of step m
%                              REACT(i,2,m) = force in Y direction at end of step m
%                              REACT(i,3,m) = moment about Z direction at end of step m
%       ELE_FOR(i,1:6,m) ==  element i's calculated 6 element end forces for all m steps
%                            of the analysis.  If this is a new analysis, ELE_FOR should
%                            be set equal to [].  If this is a continuing analysis,
%                            ELE_FOR should be first set equal to ele_for,
%                            i.e. ELE_FOR = ele_for;
%                            Note: All values reference the element's local
%                                  coordinate system.
%                            Example: ELE_FOR(17,4,5) is element 17's axial or X-
%                                     force at the end node of the element at the
%                                     end of the fifth load increment.
%                              ELE_FOR(i,1,m)  = x-force at start node at end of step m
%                              ELE_FOR(i,2,m)  = y-force at start node at end of step m
%                              ELE_FOR(i,3,m)  = z-moment at start node at end of step m
%                              ELE_FOR(i,4,m)  = x-force at end node at end of step m
%                              ELE_FOR(i,5,m)  = y-force at end node at end of step m
%                              ELE_FOR(i,6,m)  = z-moment at end node at end of step m
%       ELE_YLD(i,1:2,m) ==  applied load ratio at which element i's start and/or end node
%                            formed a plastic hinge for all m steps of the analysis.  If
%                            this is a new analysis, ELE_YLD should be set equal to [].
%                            If this is a continuing analysis, ELE_YLD should be first set
%                            equal to ele_yld, i.e. ELE_YLD = ele_yld;
%                            Be sure to insert 0.0's if at the end of a load step a 
%                            plastic hinge has not formed at the start or end of
%                            the element.
%                            Example: ELE_YLD(13,1,14) is the applied load ratio for which
%                                     a hinge formed at the start node of element 13 (as
%                                     of the end of step 14).  The value may be 0.0 which
%                                     indicates no hinge or it may be some positive real
%                                     number indicating the applied load ratio when the
%                                     hinge most recently formed. Note: if the element end
%                                     yields and then unloads elastically, then this value
%                                     should be set back to 0.0.
%                              ELE_YLD(i,1,m) = applied load ratio when a plastic hinge
%                                              first formed at the start node of element i
%                              ELE_YLD(i,2,m) = applied load ratio when a plastic hinge
%                                              first formed at the end node of element i
%       AFLAG            ==  flag to indicate if a successful analysis has or has not been
%                            completed
%                              AFLAG = 1     Successful
%                              AFLAG = 0     Unstable Structure
%                              AFLAG = -1    Analysis Halted: Limit Load Reached
%                              AFLAG = 2     Analysis Halted: Extreme Deflections
%                              AFLAG = 4     Analysis Halted: Unload/Reload Problems
%                              AFLAG = inf   No analysis code available
%       APRATIOS(i,1)    ==  applied load ratio at end of load step i
%                            Note that APRATIOS is a COLUMN vector.
%                            size(APRATIOS,1) = total number of load steps completed.
%                            If this is a new analysis, APRATIOS should be
%                            set equal to [].  If this is a continuing analysis, APRATIOS
%                            should be first set equal to apratios,
%                            i.e. APRATIOS = apratios;
%       LIMIT_STATE      ==  flag indicating post-limit state analysis
%                              LIMIT_STATE = 0  limit state not reached, system is loading
%                              LIMIT_STATE = 1  limit state reached or exceeded, system is
%                                               unloading
%
%       Version 1.0/Student's Initials/Date of Modification
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Start by defining all output arrays
%
	if restart == 1  |  isempty(apratios)
		DEFL=[]; REACT=[]; ELE_FOR=[]; ELE_YLD=[]; APRATIOS=[]; LIMIT_STATE=0;
	else
		DEFL = defl;
		REACT = react;
		ELE_FOR = ele_for;
		ELE_YLD = ele_yld;
		APRATIOS = apratios;
		LIMIT_STATE = limit_state;
	end
%
	AFLAG = inf;
%

%  STUDENT NOTE:
%     In order for this routine to become fully active AFLAG
%     must be changed.
%
%
%  Student's code starts here...
% Instantiate an object of the Analysis class
    Analysis= CMDT_Analysis_2d2in(nnodes, coord, fixity, concen, ...
                        nele, ends, A, Ayy, Izz, E, v, truss, numsteps, ...
                        ratio_req, stop_ratio,Fy,Zzz);

% % PLot Error and load norm indices
%     figure ;
%     Analysis.PlotNorms(numsteps,stop_ratio,ratio_req);

    % Extract the matrices to be returned to Mastan2
    [DEFL, REACT, ELE_FOR, AFLAG, APRATIOS, LIMIT_STATE,ELE_YLD] = Analysis.GetMastan2Returns();
    ELE_YLD
% save('CEEProgrammingAssignmentTestFrame.mat');
%
