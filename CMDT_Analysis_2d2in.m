% Cesar Y. Marco & David Tse (CMDT)
% CEE 282: Programming Project
% 2d2in Analysis Class
% 03/18/16

classdef CMDT_Analysis_2d2in < RC_Analysis_2d1el
% This is the child class of RC_Analysis_2d1el.m.
% All the protected & public data properties of parent class will be 
% inherited to this class 

% Analysis class for 2nd order analysis of a 2-dimensional structure
        properties (Access = public)
%         Miscellaneous
        Apratio
        tau_min
        
%         Tolerances
        tau_min_tol
        radtolerance
        tolerance
%         % Resultant: Free DOF
        R_free
%         
%         % New MASTAN2 returns
        APRATIOS
        LIMIT_STATE
%         
%         % Incremental version of MASTAN outputs
        DEFLstep
        REACTstep
        ELE_FORstep
        Kffstep
        ELE_YLD
        ELE_YLDstep
%         % Used for calculating error or Norm indices
        Error
        Load_Norm
        Energy_Norm
        
%         % Ratio required is often used, so create a protected property
        ratio_req

    end
    %% Public methods
    methods (Access = public)
        % Constructor inherits properties from the 2d1el analysis class
        % Add 2d2el specific object parameters, including:
        % 1) numsteps
        % 2) ratio_req
        % 3) stop_ratio
        % 4) restart
        function self = CMDT_Analysis_2d2in(nnodes, coord, fixity, concen, ...
                        nele, ends, A, Ayy, Izz, E, v, truss, numsteps, ...
                        ratio_req, stop_ratio,Fy,Zzz)
           self = self@RC_Analysis_2d1el(nnodes, coord, fixity, concen, ...
               nele, ends, A, Ayy, Izz, E, v, truss,Fy,Zzz);

           % Considers whether yielding has occurred
           self.tolerance = 0.01;
           % Consider with error in the regula falsi
           self.tau_min_tol = 1e-10;
           % Consider for radial return
           self.radtolerance = 0.03;
           
           self.ratio_req = ratio_req;
           self.Apratio = 0;
           % Initialize output variables and other necessary variables
           self.InitializeOutputVariables(); % CMDT_Analysis_2d2el
           self.CreateLoadVectors(); % RC_Analysis_2d1el
           
           %Computes initial stiffness matrix with kg = 0 and ke = #
           self.CreateStiffnessMatrix(1) % RC_Analysis_2d1el
           
           % Computes DEFLstep and REACTstep
           self.ComputeDisplacementsReactions(self.Pf*ratio_req);% CMDT_Analysis_2d2el
           
           % Run the incremental single step analysis with no error
           % iterations
           self.RunAnalysis(numsteps,stop_ratio,ratio_req); % CMDT_Analysis_2d2el
        end
        
        %% Run the single step incremental analysis with no error iterations
        function RunAnalysis(self,numsteps,stop_ratio,ratio_req)
            % Run a for loop to start at 1 & end at a user-specified point
            % This for loop the minimum of the number of steps and stop
            % ratio divided by the ratio req.
%             endpoint = (min(numsteps,stop_ratio/ratio_req));
            i = 1;
            while i <= numsteps && self.Apratio <= stop_ratio
                % for each element calculate the local geometric stiffness
                % matrix to be added with ke during each load step
                for j = 1:self.nele
                    self.elements(j).ComputeLocalGeometricStiffnessMatrix()
                    % CMDT_Element_2d2el
                end
%                 
                % Recreate the stiffness matrix based on the new
                % geometry, which would affect only the kg and the
                % transformation matrix. The local ke would not change
                % because this value is only calculated once
                self.CreateStiffnessMatrix(i)
                
                % Run the analysis only if the Kff matri is
                % well-conditioned and K is positive definite according to
                % the results from cholesky decomposition
                if (self.AFLAG == 1) && (self.LIMIT_STATE == 0)
                    % Calculate the global DEFL and REACT to be used for
                    % obtaining the element forces via natural deformation
                    self.ComputeDisplacementsReactions(self.Pf*ratio_req); % RC_Analysis_2d1el
                    % Calculate the theta_an,theta_bn and un used to
                    % calculate the incremental element forces
                    self.RecoverElementForces(i); % CMDT_Analysis_2d1in
                    %
                    %
                    %
                    %% RETURN METHOD
                    self.tau_min = 1;
                    for elenum = 1:self.nele
                        
                        % Obtain the incremental element force to determine
                        % the P-M values as well as the incremental dP and
                        % dM values
                        dF = self.elements(elenum).GetdFlocal();
                        End_i_dP = dF(1);
                        End_i_dM = dF(3);
                        End_j_dP = dF(4);
                        End_j_dM = dF(6);
                        if i > 1
                            F = self.ELE_FOR(elenum,:,i-1);
                            End_i_P = F(1);
                            End_i_M = F(3);
                            End_j_P = F(4);
                            End_j_M = F(6);
                        else
                            End_i_P = 0;
                            End_i_M = 0;
                            End_j_P = 0;
                            End_j_M = 0;
                        end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% i-th end calculations when the PM location is below the YS and is going to
% breach the YS
                        % Compute the prior step and post step to determine
                        % which case the step falls in
                        Yield_Surface_End_i_prior_step = ...
                            self.elements(elenum).ComputeYieldSurface(End_i_P,End_i_M);
                        Yield_Surface_End_i_post_step = ...
                            self.elements(elenum).ComputeYieldSurface(...
                            End_i_P + End_i_dP,...
                            End_i_M + End_i_dM);
                        self.Calc_Tau_min (1,Yield_Surface_End_i_prior_step,...
                            Yield_Surface_End_i_post_step, ...
                            self.tolerance, 1.0 + self.tolerance,...
                            End_i_P, End_i_M, End_i_dP, End_i_dM, elenum, i);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% j-th end calculations when the PM location is below the YS and is going to
% breach the YS
                        Yield_Surface_End_j_prior_step = ...
                            self.elements(elenum).ComputeYieldSurface(End_j_P,End_j_M);
                        Yield_Surface_End_j_post_step = ...
                            self.elements(elenum).ComputeYieldSurface(...
                            End_j_P + End_j_dP,...
                            End_j_M + End_j_dM);
                        self.Calc_Tau_min (2,Yield_Surface_End_j_prior_step,...
                            Yield_Surface_End_j_post_step,...
                            self.tolerance, 1.0 + self.tolerance, End_j_P, ...
                                End_j_M, End_j_dP, End_j_dM, elenum, i);       
                    end
                    %Update the Applied load ratio
                    if i > 1
                        if self.APRATIOS(i-1)+self.ratio_req > ...
                                self.APRATIOS(i-1)+self.tau_min*self.ratio_req
                           self.APRATIOS(i)=self.APRATIOS(i-1)+self.tau_min*self.ratio_req;
                        else
                            self.APRATIOS(i)=self.APRATIOS(i-1)+self.ratio_req;   
                        end
                    else
                        self.APRATIOS(i)=self.ratio_req;
                    end
                    self.Apratio = self.APRATIOS(i);
                    % algorithm to format what hinges yield 
                     ele_yield_step = find(self.ELE_YLDstep > self.APRATIOS(i));
                    self.ELE_YLDstep(ele_yield_step) = 0;        
                    % Recalculate the displacement, reactions, and element
                    % forces by including the tau_min coefficient
                    % calculated earlier
                    if self.tau_min < 1
                        self.ComputeDisplacementsReactions(self.Pf * ratio_req * self.tau_min);
                        self.RecoverElementForces(i);
                    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOADING AN ELEMENT ALREADY ON THE YS TO ANOTHER POINT ON THE YS
                    for elenum = 1:self.nele
                        dF = self.elements(elenum).GetdFlocal();
                        End_i_dP = dF(1);
                        End_i_dM = dF(3);
                        End_j_dP = dF(4);
                        End_j_dM = dF(6);
                        
                        if i > 1
                            F = self.ELE_FOR(elenum,:,i-1);
                            End_i_P = F(1);
                            End_i_M = F(3);
                            End_j_P = F(4);
                            End_j_M = F(6);
                        else
                            End_i_P = 0;
                            End_i_M = 0;
                            End_j_P = 0;
                            End_j_M = 0;
                        end
                        Yield_Surface_End_i_post_step = ...
                            self.elements(elenum).ComputeYieldSurface(...
                            End_i_P + End_i_dP,...
                            End_i_M + End_i_dM);
                        Yield_Surface_End_j_post_step = ...
                            self.elements(elenum).ComputeYieldSurface(...
                            End_j_P + End_j_dP,...
                            End_j_M + End_j_dM);
                        if Yield_Surface_End_i_post_step > (1 + self.radtolerance)
                            self.RadialReturn(End_i_P, End_i_M, End_i_dP, ...
                                End_i_dM, elenum, [1;3]);
                        end
                        if Yield_Surface_End_j_post_step > (1 + self.radtolerance)
                            self.RadialReturn(End_j_P, End_j_M, End_j_dP, ...
                                End_j_dM, elenum, [4;6]);
                        end
                    end
                % for steps 2 onwards, take a cumulative sum approach
                % for calculating the global deflections and reactions
                if(i > 1)
                    self.DEFL(:,:,i) = cat(3,self.DEFLstep+self.DEFL(:,:,i-1));
                    self.REACT(:,:,i) = cat(3,self.REACTstep+self.REACT(:,:,i-1));
                    self.ELE_FOR(:,:,i) = cat(3,self.ELE_FORstep);
                    self.ELE_YLD(:,:,i) = cat(3,self.ELE_YLDstep);
                    % for step 1 just add the increment to the initial zero
                    % 3D matrix
                else
                    self.ELE_YLD(:,:,1) = self.ELE_YLD(:,:,1) + self.ELE_YLDstep;
                    self.DEFL(:,:,1) =  self.DEFL(:,:,1) + self.DEFLstep;
                    self.REACT(:,:,1) = self.REACT(:,:,1) + self.REACTstep;
                    self.ELE_FOR(:,:,1) = self.ELE_FOR(:,:,1) + self.ELE_FORstep;
                end

            % Update the node coordinates
            self.UpdateGeometry(); % CMDT_Analysis_2d2el
                        
            % Computations used for error and norm
            self.ComputeResultant(); % CMDT_Analysis_2d2el
            self.ComputeError(i); % CMDT_Analysis_2d2el
            self.ComputeNorms(i); % CMDT_Analysis_2d2el
            i = i + 1;
                
                end
            if self.LIMIT_STATE == 1
                % break the for loop if the above if condition is not
                % met, this means that either a limit is reached or the
                % stiffness matrix is ill-conditioned
                break
            end
            end    
        end
        %% Calulate the tau_min used to scale increments that breach the yield surface
        function Calc_Tau_min (self,hinge,YSprestep,YSpoststep, YStolerance, ...
                max_tolerance, P, M, dP, dM, elementnum, stepnum)
% hinge:
    % 1 = ith-end
    % 2 = jth-end
% YSprestep: Prior step PM location of the Yield Surface
% YSpoststep: Post step PM location of the Yield Surface
% YStolerance: tolerance for the YS of 1.0
% max_tolerance: max tolerable YS
% P: axial force
% M: Moment
% dP: incremental axial force
% dM: incremental moment
% elenum: element number
% stepnum: step number


% Case A: When the YS post step has not breached the YS
                if YSpoststep < (1 - YStolerance)
                    self.ELE_YLDstep(elementnum,hinge) = 0;   
% Case B: When the YS poststep is between 1.0-YS_tolerance and 1.0
                elseif YSpoststep >= (1 - YStolerance) && YSpoststep <= 1
                    if self.ELE_YLDstep(elementnum,hinge) == 0;             
                      self.ELE_YLDstep(elementnum,hinge) = self.APRATIOS(stepnum-1)...
                           + self.ratio_req;      
                    end
% Case C: When the YS poststep is between 1.0 and 1.0+YS_tolerance
                elseif YSpoststep > 1  && YSpoststep <= (1 + YStolerance)
                    if YSprestep < (1 - YStolerance) && self.ELE_YLDstep(elementnum,hinge)==0;             
                        tau_update = self.ComputeReturn(dP,dM,...
                         P,M,self.elements(elementnum),1.0);
                      self.ELE_YLDstep(elementnum,hinge) = self.APRATIOS(stepnum-1)...
                           + tau_update*self.ratio_req;   
                       if tau_update < self.tau_min
                        self.tau_min = tau_update;
                       end
                    end
% Case D: When the YS poststep is between 1.0+YS_tolerance and the max
% specified tolerance
                elseif YSpoststep > (1 + YStolerance) && YSpoststep <= max_tolerance
                    if YSprestep < (1 - YStolerance) && self.ELE_YLDstep(elementnum,hinge)==0;             
                    	tau_update = self.ComputeReturn(dP,dM,...
                         P,M,self.elements(elementnum),1.0);
                        if self.ELE_YLDstep(elementnum,hinge) == 0
                            self.ELE_YLDstep(elementnum,hinge) = self.APRATIOS(stepnum-1)...
                                + tau_update*self.ratio_req;  
                        end
                        if tau_update < self.tau_min
                        self.tau_min = tau_update;
                        end
                    end

% Case E: When the YS poststep is lower than the max_tolerance specified
                elseif YSpoststep >= max_tolerance
                    if YSprestep >= (1 - YStolerance)
                    tau_update = self.ComputeReturn(dP,dM,...
                        P,M,self.elements(elementnum),max_tolerance);
                    elseif YSprestep < (1 - YStolerance)            
                        tau_update = self.ComputeReturn(dP,dM,...
                         P,M,self.elements(elementnum),1.0);
                        if self.ELE_YLDstep(elementnum,hinge) == 0
                            self.ELE_YLDstep(elementnum,hinge) = self.APRATIOS(stepnum-1)...
                            + tau_update*self.ratio_req;   
                        end
                    end
                    if tau_update < self.tau_min
                        self.tau_min = tau_update;
                    end
                end
        end
        %% Radial return method for when a hinge is already yielding and strays away from the YS
        function RadialReturn(self, P, M, dP, dM, elementnum, DOF)
% P: axial force
% M: moment
% dP: incremental axial force
% dM: incremental moment
% elenum: element number
% DOF: Degrees of freedom
            tau_update = self.ComputeReturn(0 - (P+dP), 0-(M+dM), P+dP , M+dM , ...
                self.elements(elementnum) , 1.0);
            self.ELE_FORstep(elementnum , DOF(1)) = P+tau_update * (0-(P+dP));
            self.ELE_FORstep(elementnum , DOF(2)) = M+tau_update * (0-(M+dM));
        end
        
        %% Plots the load norm indices & energy norm indices with the Applied Load Ratio
        function PlotNorms(self)
            % x-axis: Applied load ratio
            
            % y-axis: energy or load norm index values
            plot(self.APRATIOS,self.Load_Norm,self.APRATIOS,self.Energy_Norm);
            
            % Plotting modifications
            xlabel('Applied Load Ratio');
            ylabel('Error Index');
            title('Error indices vs. Applied load ratio');
            legend('Load Norm Error Index','Energy Norm Error Index');
        end
        %% Get Mastan2 Returns
        %  Returns the matrices that need to be returned to Mastan2
        function [DEFL, REACT, ELE_FOR, AFLAG, APRATIOS, LIMIT_STATE,ELE_YLD] = ...
                  GetMastan2Returns(self)
            DEFL = self.DEFL;
            REACT = self.REACT;
            ELE_FOR = self.ELE_FOR;
            AFLAG = self.AFLAG;
            APRATIOS = self.APRATIOS';
            LIMIT_STATE = self.LIMIT_STATE;
            ELE_YLD = self.ELE_YLD;
        end
    end
    %% Protected Methods
    methods (Access = protected)
        %% Create Nodes
        % Create the nnodes x 1 vector of 2d2el node objects representing all 
        % the nodes in the structure
        function CreateNodes(self)
            for i = 1:self.nnodes
                self.nodes = [self.nodes; CMDT_Node_2d2in(i, self.coord_t(:,i))];
            end
        end
        %% Create Elements
        % Create the nele x 1 vector of 2d2el element objects representing all 
        % the elements in the structure
        % Arguments are all matrices received from Mastan2. Refer to 
        % comments in ud_2d1el.m for details.
        function CreateElements(self, A, Ayy, Izz, E, v, Fy, Zzz)
            for i = 1:self.nele
                % Create an Element object and append it to the "elements" vector
                self.elements = [self.elements; CMDT_Element_2d2in(...
                    self.nodes(self.ends(i, 1:2)), A(i), Ayy(i), Izz(i),...
                    E(i), v(i), self.truss,Fy(i),Zzz(i))];
            end
        end
        %% Compute Displacements Reactions
        %  Compute the displacements and reactions and format them to return to Mastan2
        function ComputeDisplacementsReactions(self,dP) % RC_Analysis_2d1el
            % Initialize
            self.DEFLstep = zeros(3,self.nnodes);
            self.REACTstep = zeros(3,self.nnodes);

            % Compute the displacements, dD
            self.delf = self.Kff \ (dP - self.Kfn*self.deln);
            
            % Format the computed displacements using linear indexing of 
            % the "DEFL" matrix
            self.DEFLstep(self.dof_free) = self.delf;
            self.DEFLstep(self.dof_disp) = self.deln;
            
            % Columns: DOF, Rows: Nodes
            self.DEFLstep = self.DEFLstep';
            
            % Compute the reactions, accounting for loads applied directly 
            % on the supports
            self.Ps = self.Ksf*self.delf + self.Ksn*self.deln - ...
                      self.Psupp * self.ratio_req ;
            self.Pn = self.Knf*self.delf + self.Knn*self.deln;
            
            % Format the computed reactions using linear indexing of the 
            % "REACT" matrix
            self.REACTstep(self.dof_supp) = self.Ps;
            self.REACTstep(self.dof_disp) = self.Pn;
            self.REACTstep = self.REACTstep';
        end

        %% Recover element forces based on the natural deformation approach
        function RecoverElementForces(self,step)
            % Columns: Nodes, Rows: DOF
            DEFL_t = self.DEFLstep';
            for i = 1:self.nele
                % Obtain the displacements at the degrees of freedom 
                % corresponding to element i using linear indexing of the 
                % "DEFL_t" matrix
                self.elements(i).ComputeForces(...
                                 DEFL_t(self.elements(i).GetElementDOF()));
                if step >1
                    self.ELE_FORstep(i,:) = self.ELE_FOR(i,:,step-1) + ...
                                        self.elements(i).GetdFlocal()';
                else
                    self.ELE_FORstep(i,:) = self.ELE_FOR(i,:,1) + ...
                                        self.elements(i).GetdFlocal()';
                end
                self.elements(i).UpdateFLocal(self.ELE_FORstep,i);   
            end
        end
        %% Global Coordinates, updates after the load step
        % Purpose: Update the geometry to calculate the incremental output
        % variables, such as the incremental internal force
        function UpdateGeometry(self)
            % Columns: DOF, Rows: Nodes
            defl_t = self.DEFLstep';
            for i = 1:self.nnodes
                nodeDOF = self.nodes(i).GetNodeDOF();
                % Obtain only the displacements in the x and y direction to
                % update the geometry
                nodeDOF = nodeDOF(1:2);
                
                self.nodes(i).UpdateCoordinates(self.nodes(i).GetNodeCoord()...
                    + defl_t(nodeDOF));
            end
        end
        %% Computes the Resultant
        % Purpose: The resultant is used to calculate the error where 
        % {E}={P}-{R}
        function ComputeResultant(self)
            R=zeros(self.num_dof_total,1);
    
            %loop over all elements, pulling each elements DOF and updating their 
            %transformation matrix. Use DOF to ensure that the
            %corresponding global force of each element is added to the
            %corresponding entry in the resultant vector.
            for i = 1:self.nele
               %pull current elements dof
               element_dof = self.elements(i).GetElementDOF();
              self.elements(i).UpdateTransformationMatrix(); 
               gamma = self.elements(i).GetTransformationMatrix();
               %update element forces to global coordinates
               F_global = gamma'*self.ELE_FORstep(i,:)'; 
               
               for j=1:length(element_dof)
                   % add the global forces of the element to Resultant
                   % based on the corresping dof
                   R(element_dof(j))=R(element_dof(j))+F_global(j); 
               end 
            end
            self.R_free=R(self.dof_free);
        end

        %% Compute Error
        % Computes the incremental error based on the applied load and
        % resultant
        function ComputeError(self,step)
            self.Error = self.Pf*self.APRATIOS(step) - self.R_free;
        end
        %% Compute load norm indices and the error indices
        % Purpose: Compute norms to be plotted against the applied load
        % ratio
        function ComputeNorms(self,i)
            self.Load_Norm = [self.Load_Norm ; norm(self.Error) / ...
                              norm(self.Pf*self.APRATIOS(i))];
            self.Energy_Norm = [self.Energy_Norm;(abs(self.Error') * ...
                                abs(self.delf))/(abs(self.Pf'*self.APRATIOS(i))...
                                *abs(self.delf))];
        end
        %% Check Kff Matrix
        % Purpose: Determines whether the incremental SS analysis should stop 
        % prematurely or reaches the user-specified end points based on
        % matrix conditioning and cholesky decomposition
        
        % AFLAG = 1     Successful
        % AFLAG = 0     Unstable Structure
        % AFLAG = -1    Analysis Halted: Limit Load Reached
        % LIMIT_STATE = 0  limit state not reached, system is loading
        % LIMIT_STATE = 1  limit state reached or exceeded, system is
        %                  unloading
        function CheckKffMatrix(self)
            [~,p] = chol(self.Kff);
            if self.Kffstep == 1 && ((condest(self.Kff) ...
               > self.Kff_condition_threshold) || (p > 0))
                self.LIMIT_STATE = 0;                
                self.AFLAG = 0;
            elseif (p > 0) % Limit Point based on Cholesky decomposition
                self.LIMIT_STATE = 1;
                self.AFLAG = -1;
            else
                self.LIMIT_STATE = 0;
                self.AFLAG = 1;
            end
            self.Kffstep = self.Kffstep + 1;
        end
        %% Initialize variables to be able to run the analysis
        % 2D initialization does not include the number of steps
        function InitializeOutputVariables(self)
           self.Kffstep = 1;
           self.DEFL = zeros(self.nnodes,self.num_dof_node,1);
           self.REACT = zeros(self.nnodes,self.num_dof_node,1);
           self.ELE_FOR = zeros(self.nele,self.num_dof_node*2,1);
           self.ELE_FORstep = zeros(self.nele, self.num_dof_node*2);
           self.ELE_YLD = zeros(self.nele,2,1);
           self.ELE_YLDstep = zeros(self.nele, 2);
        end
        %% Reguli falsi method
        % Return an attempted load increment to below the yield surface
        function [tau_regula] = ComputeReturn(self,dP,dM,P,M,element,max_limit)
            % Pass in YS_current_step and YS_next_step for one end at a time, 
            % i.e. a 2x1 vector. 
            % one element and one hinge at a time
            % Run if statement to determine the state of the element, it
            % could be 1 out of 3 cases, set tau = 0, and only update it
            % with tau_new if tau_new is lower in value than tau       
            
            % Initialize
            tau_lower = 0;
            tau_upper = 1;
            error = 10;
            % Error reduction based on while loop
            while error > self.tau_min_tol
                phi_lower = element.ComputeYieldSurface(P+tau_lower*dP,...
                            M+tau_lower*dM);
                phi_upper = element.ComputeYieldSurface(P+tau_upper*dP,...
                            M+tau_upper*dM);
            	tau_regula = tau_upper - ((phi_upper-max_limit)*(tau_lower-tau_upper))/...
                             (phi_lower - phi_upper);
                phi_regula = element.ComputeYieldSurface(P+tau_regula*dP ,...
                             M+tau_regula*dM);
                if sign(phi_regula - max_limit) == sign(phi_lower - max_limit)
                    tau_lower = tau_regula;
                else 
                    tau_upper = tau_regula;   
                end
                error = abs(max_limit - phi_regula);
            end
        end
        %% Create Stiffness Matrix
        %  Create the global stiffness matrix for the structure and store 
        %  it in sparse format modified to include the gradient and yield 
        % surface calculations
        function CreateStiffnessMatrix(self,step)
            
            % Initialize the vectors that will be store the coordinates and 
            % values of the non-zero elements in K
            K_row = [];
            K_col = [];
            K_data = [];
            
            % Loop over all elements and append the contribution of each 
            % element's global stiffness matrix to the "K_row", "K_col", 
            % and "K_data" vectors
            for i = 1:self.nele
                if step > 1
                    F = self.ELE_FOR(i,:,step-1);
                else
                    F = self.ELE_FOR(i,:,1);
                end
                YS_i = self.elements(i).ComputeYieldSurface(F(1),F(3));
                YS_j = self.elements(i).ComputeYieldSurface(F(4),F(6));
                self.elements(i).ComputeGradient([F(1);F(4)],[F(3);F(6)],...
                    [YS_i;YS_j],self.tolerance);
                self.elements(i).ComputeGlobalStiffnessMatrix();
                [row, col, data] = find(self.elements(i).GetKGlobal());
                
                element_dof = self.elements(i).GetElementDOF();
                K_row = [K_row; element_dof(row)];
                K_col = [K_col; element_dof(col)];
                K_data = [K_data; data];
            end
            % Convert the "K_row", "K_col", and "K_data" vectors to a sparse matrix
            self.K = sparse(K_row, K_col, K_data, self.num_dof_total, self.num_dof_total);            
            self.ComputeStiffnessSubMatrices();
            self.CheckKffMatrix();
        end     
    end
    
end
