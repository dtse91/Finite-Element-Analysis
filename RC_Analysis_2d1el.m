classdef RC_Analysis_2d1el < handle

% Analysis class for 1st order analysis of a 2-dimensional structure
    
    properties (Access = public)
        % Matrices received from Mastan2. Refer to comments in ud_2d1el.m for details.
        nnodes
        nele
        ends
        truss
        
        % Transposes of the corresponding matrices received from Mastan2. Refer to comments in ud_2d1el.m for 
        % details.
        coord_t
        fixity_t
        concen_t
        
        % Total number of degrees of freedom in the structure
        num_dof_total
        
        % Total number of degrees of freedom that are free, have specified displacements, and are supported
        num_dof_free
        num_dof_disp
        num_dof_supp
        
        % Vectors of the free, displaced, and support degree of freedom numbers
        dof_free
        dof_disp
        dof_supp
        
        % nnodes x 1 vector of node objects representing all the nodes in the structure
        nodes
        
        % nele x 1 vector of element objects representing all the elements in the structure
        elements
        
        % Global stiffness matrix for the structure (sparse)
        K
        
        % Sub-matrices of K (sparse)
        Kff
        Kfn
        Knf
        Knn
        Ksf
        Ksn
        
        % Sub-vectors of the force vector
        Pf
        Pn
        Ps
        
        % Vector of forces applied directly on the supports
        Psupp
        
        % Sub-vectors of the displacement vector
        delf
        deln
        
        % Matrices to be returned to Mastan2. Refer to comments in ud_2d1el.m for details.
        DEFL
        REACT
        ELE_FOR
        AFLAG
    end
    
    properties (Constant)
        % Number of degrees of freedom per node
        num_dof_node = 3
        
        % Largest acceptable condition number for the Kff matrix, beyond which, the structure is considered
        % unstable
        Kff_condition_threshold = 1e12
    end
    
    methods (Access = public)
        %% Constructor
        %    Arguments are all matrices received from Mastan2. Refer to comments in ud_2d1el.m for details.
        function self = RC_Analysis_2d1el(nnodes, coord, fixity, concen, nele, ends, A, Ayy, Izz, E, v, truss, Fy, Zzz)
            self.num_dof_total = nnodes*self.num_dof_node;
            
            self.nnodes = nnodes;
            self.coord_t = coord';
            self.fixity_t = fixity';
            self.concen_t = concen';
            self.nele = nele;
            self.ends = ends;
            self.truss = truss;
            
            self.CreateNodes();
            self.CreateElements(A, Ayy, Izz, E, v, Fy, Zzz);
            self.ClassifyDOF();
        end
        
        %% Run Analysis
        %  Run 1st order analysis
        function RunAnalysis(self)
            self.InitializeOutputVariables();
            
            % Create the structure's global stiffness matrix, check the condition number of the Kff matrix and
            % set AFLAG accordingly
            self.CreateStiffnessMatrix();
            
            % Run the analysis only if the structure is stable, i.e. the Kff matrix is well conditioned
            if self.AFLAG
               self.CreateLoadVectors();
               self.ComputeDisplacementsReactions();
               self.RecoverElementForces();
            end
        end
        
        %% Get Mastan2 Returns
        %  Returns the matrices that need to be returned to Mastan2
        function [DEFL, REACT, ELE_FOR, AFLAG] = GetMastan2Returns(self)
            DEFL = self.DEFL;
            REACT = self.REACT;
            ELE_FOR = self.ELE_FOR;
            AFLAG = self.AFLAG;
        end
    end
    
    methods (Access = protected)
        %% Create Nodes
        %  Create the nnodes x 1 vector of node objects representing all the nodes in the structure
        function CreateNodes(self)
            for i = 1:self.nnodes
                % Create a Node object and append it to the "nodes" vector
                self.nodes = [self.nodes; RC_Node_2d1el(i, self.coord_t(:,i))];
            end
        end

        %% Create Elements
        %  Create the nele x 1 vector of element objects representing all the elements in the structure
        %    Arguments are all matrices received from Mastan2. Refer to comments in ud_2d1el.m for details.
        function CreateElements(self, A, Ayy, Izz, E, v)
            for i = 1:self.nele
                
                % Create an Element object and append it to the "elements" vector
                self.elements = [self.elements; RC_Element_2d1el(self.nodes(self.ends(i, 1:2)), A(i), ...
                                    Ayy(i), Izz(i), E(i), v(i), self.truss)];
            end
        end
        
        %% Classify DOF
        %  Create vectors of degrees of freedom that are free, have specified displacements, and correspond to
        %  supports. Also compute the total number of degrees of freedom in each category.
        function ClassifyDOF(self)
            
            % Use linear indexing of the "fixity_t" matrix to obtain each vector of degrees of freedom
            self.dof_free = find(isnan(self.fixity_t));
            self.num_dof_free = length(self.dof_free);
            
            self.dof_disp = find(self.fixity_t ~= 0 & ~isnan(self.fixity_t));
            self.num_dof_disp = length(self.dof_disp);
            
            self.dof_supp = find(self.fixity_t == 0);
            self.num_dof_supp = length(self.dof_supp);
        end
        
        %% Initialize Output Variables
        %  Initialize the matrices to be returned to Mastan2 with zeros
        function InitializeOutputVariables(self)
            self.DEFL = zeros(self.num_dof_node, self.nnodes);
            self.REACT = zeros(self.num_dof_node, self.nnodes);
            self.ELE_FOR = zeros(self.nele, self.num_dof_node*2);
        end
        
        %% Create Stiffness Matrix
        %  Create the global stiffness matrix for the structure and store it in sparse format
        function CreateStiffnessMatrix(self)
            
            % Initialize the vectors that will be store the coordinates and values of the non-zero elements
            % in K
            K_row = [];
            K_col = [];
            K_data = [];
            
            % Loop over all elements and append the contribution of each element's global stiffness matrix to 
            % the "K_row", "K_col", and "K_data" vectors
            for i = 1:self.nele
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
        
        %% Compute Stiffness Sub Matrices
        %  Compute the partitions of the K matrix
        function ComputeStiffnessSubMatrices(self)
            self.Kff = self.K(self.dof_free, self.dof_free);
            self.Kfn = self.K(self.dof_free, self.dof_disp);
            self.Knf = self.K(self.dof_disp, self.dof_free);
            self.Knn = self.K(self.dof_disp, self.dof_disp);
            self.Ksf = self.K(self.dof_supp, self.dof_free);
            self.Ksn = self.K(self.dof_supp, self.dof_disp);
        end
        
        %% Check Kff Matrix
        %  If the condition number of the Kff matrix is greater than the allowable threshold, set AFLAG to 0
        %  to indicate the structure is unstable, else set it to 1 to continue running the analysis
        function CheckKffMatrix(self)
            if condest(self.Kff) > self.Kff_condition_threshold
                self.AFLAG = 0;
            else
                self.AFLAG = 1;
            end
        end
        
        %% Create Load Vectors
        %  Create the applied load vectors
        function CreateLoadVectors(self)
            
            % Compute vector of concentrated loads applied at the free and support degrees of freedom using
            % linear indexing of the "concen_t" matrix
            self.Pf = self.concen_t(self.dof_free);
            self.Psupp = self.concen_t(self.dof_supp);
            
            % Compute vector of specified displacements using linear indexing of the "fixity_t" matrix
            self.deln = self.fixity_t(self.dof_disp);
        end
        
        %% Compute Displacements Reactions
        %  Compute the displacements and reactions and format them to return to Mastan2
        function ComputeDisplacementsReactions(self)
            
            % Compute the displacements
            self.delf = self.Kff \ (self.Pf - self.Kfn*self.deln);
            
            % Compute the reactions, accounting for loads applied directly on the supports
            self.Ps = self.Ksf*self.delf + self.Ksn*self.deln - self.Psupp;
            self.Pn = self.Knf*self.delf + self.Knn*self.deln;
            
            % Format the computed displacements using linear indexing of the "DEFL" matrix
            self.DEFL(self.dof_free) = self.delf;
            self.DEFL(self.dof_disp) = self.deln;
            self.DEFL = self.DEFL';
            
            % Format the computed reactions using linear indexing of the "REACT" matrix
            self.REACT(self.dof_supp) = self.Ps;
            self.REACT(self.dof_disp) = self.Pn;
            self.REACT = self.REACT';
        end
        
        %% Recover Element Forces
        %  Recover the local element forces and format them to return to Mastan2
        function RecoverElementForces(self)
            DEFL_t = self.DEFL';
            for i = 1:self.nele
                
                % Obtain the displacements at the degrees of freedom corresponding to element i using linear
                % indexing of the "DEFL_t" matrix
                self.elements(i).ComputeForces(DEFL_t(self.elements(i).GetElementDOF()));
                self.ELE_FOR(i,:) = self.elements(i).GetFLocal();
            end
        end
    end
end