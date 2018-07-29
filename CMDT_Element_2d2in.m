% Cesar Y. Marco & David Tse (CMDT)
% CEE 282: Programming Project
% 2d2in Element Class
% 03/18/16

classdef CMDT_Element_2d2in < RC_Element_2d1el
% Replace XYZ by your initials and rename the file accordingly before proceeding
% This is the child class of RC_Element_2d1el.m.
% All the protected & public data properties of parent class will be inherited to this class 


% Element class for 2nd order analysis of a 2-dimensional structure
    
    % Protected properties go here
    properties (Access = public)
        % Use this space to define all the additional dataproperties if required.
        kg
        K
        dFLocal
        Fy
        Mp
        G
        YS
        kp
    end
    
    % Public methods go here
    methods (Access = public)
        % Define the constructor here and use it to call the parent class.
        function self = CMDT_Element_2d2in(element_nodes, A, Ayy, Izz, E, v, truss,Fy,Zzz)
           self = self@RC_Element_2d1el(element_nodes, A, Ayy, Izz, E, v, truss);
           self.dFLocal = zeros(6,1);
           self.f_local = zeros(6,1);
           self.kp = 0;
           self.ComputeLocalGeometricStiffnessMatrix();
           syms P M
           self.ComputeYieldStrength(Fy,A);
           self.ComputePlasticMoment(Zzz,Fy);
           self.YS = self.ComputeYieldSurface(P,M)
        end       
        %% Compute the geometric stiffness matrix in global coordinates
        function ComputeGlobalStiffnessMatrix(self)
            self.K = self.gamma' * (self.ke_local + self.kp + self.kg)*self.gamma;
        end
        %% Getter Functions
        % Allows superclasses (i.e. 2d2el Analysis) to obtain protected
        % variables from this class
        function K = GetKGlobal(self)
            K = self.K;
        end
        
        function dF = GetdFlocal(self)
            dF = self.dFLocal;
        end    
        
        function gamma = GetTransformationMatrix(self)
            gamma = self.gamma;
        end
        %% Updates the transformation matrix based on the new geometry
        % Used in ComputeResultant() in the 2d2el analysis class 
        function UpdateTransformationMatrix(self)
            axis = self.element_nodes(2).GetNodeCoord() - self.element_nodes(1).GetNodeCoord();
            % Calculate terms
            theta = cart2pol(axis(1), axis(2));
            c = cos(theta);
            s = sin(theta);
            % Assemble 3x3 gamma terms
            gamma3 = [ c,  s,  0;
                      -s,  c,  0;
                       0,  0,  1];
            zeros3 = zeros(3);
            % Assemble big gamma
            self.gamma = sparse([gamma3, zeros3;
                                 zeros3, gamma3]);
        end 
        %%  Updates the element force
        % This method is used within the the RecoverElementForces() method
        % in the 2d2el analysis class
        function UpdateFLocal(self,ELE_FOR,ele_num)
            self.f_local = ELE_FOR(ele_num,:)';
        end
        %% Element forces based on the natural deformation approach
         function ComputeForces(self, del_global)
            self.del_global = del_global;
            
            % Compute the element displacement vector in local coordinates
            self.del_local = self.gamma * self.del_global;
            
            % Compute theta_r
            theta_r = atan((self.del_local(5) - self.del_local(2))/...
                      (self.L + self.del_local(4) - self.del_local(1))); %rad
                  
            % Compute theta_an
            theta_an = self.del_local(3) - theta_r;
            
            % Compute theta_bn
            theta_bn = self.del_local(6) - theta_r;
            
            % Compute un 
            un = (self.del_local(4)-self.del_local(1)) ...
                + ((self.del_local(4)-self.del_local(1))^2 + ...
                (self.del_local(5)-self.del_local(2))^2)/(2*self.L);
            
            % Assemble del_delta_n in local coordinates
            dDn = [0;0;theta_an;un;0;theta_bn];
            
%             Compute the element force vector in local coordinates
            self.dFLocal = (self.ke_local + self.kp + self.kg) * dDn;
%             self.dFLocal = (self.ke_local + self.kp) * self.del_local;

         end
        %% Compute Local Geometric Stiffness Matrix
        %  Check whether the element is part of a truss or not, and compute its local elastic stiffness matrix
        %  accordingly. Store the computed matrix in sparse format.
        function ComputeLocalGeometricStiffnessMatrix(self)
            self.kg = sparse(6,6);
            self.kg(1,1) = self.f_local(4)/self.L;
            self.kg(1,4) = -self.kg(1,1);
            self.kg(2,2) = 1.2*self.f_local(4)/self.L;
            self.kg(2,3) = self.f_local(4)/10;
            self.kg(2,5) = -self.kg(2,2);
            self.kg(2,6) = self.kg(2,3);
            self.kg(3,3) = 2*self.L*self.f_local(4)/15;
            self.kg(3,5) = -self.kg(2,3);
            self.kg(3,6) = -self.f_local(4)*self.L/30;
            self.kg(4,4) = self.kg(1,1);
            self.kg(5,5) = self.kg(2,2);
            self.kg(5,6) = -self.kg(2,3);
            self.kg(6,6) = self.kg(3,3);
            self.kg(3,2) = self.kg(2,3);
            self.kg(4,1) = self.kg(1,4);
            self.kg(5,2) = self.kg(2,5);
            self.kg(5,3) = self.kg(3,5);
            self.kg(6,2) = self.kg(2,6);
            self.kg(6,3) = self.kg(3,6);
            self.kg(6,5) = self.kg(5,6);
        end
        %% Compute the Yield Strength
        function ComputeYieldStrength(self,fy,A)
            % Convert stress to force
            self.Fy = fy*A;
        end
        %% Compute the Plastic Moment
        function ComputePlasticMoment(self,Zx,fy)
            self.Mp = Zx*fy;
        end
        %% Compute the Yield Surface
        function YieldS = ComputeYieldSurface(self,P,M)
            YieldS = (P/self.Fy)^2 + (M/self.Mp)^2 + 3.5*(P/self.Fy)^2 * (M/self.Mp)^2;
        end
        %% Compute the gradient matrix of the yield surface
        function ComputeGradient(self,P,M,Yield,tolerance)
%             Force and Moment are 2x1 vectors, where the 1st row is the i
%             end and the 2nd row is the j end
            % Yield is a 2x1 vector where i is the 1st row and j is the 2nd
            % row, check using if statments to determine if both or one end 
            % has yielded 
            if Yield(1) >= (1 - tolerance) && Yield(2) >= (1-tolerance)
                % both i and j end yield
                self.G = zeros(6,2);
                 self.G = [ 
(2*P(1))/self.Fy^2 + (7*M(1)^2*P(1))/(self.Fy^2*self.Mp^2),                                              0;
                                                0,                                                       0;
(2*M(1))/self.Mp^2 + (7*M(1)*P(1)^2)/(self.Fy^2*self.Mp^2),                                              0;
                                                0, (2*P(2))/self.Fy^2 + (7*M(2)^2*P(2))/(self.Fy^2*self.Mp^2);
                                                0,                                                       0;
                                                0, (2*M(2))/self.Mp^2 + (7*M(2)*P(2)^2)/(self.Fy^2*self.Mp^2)]  ;
                self.ComputeKPlastic();
            elseif Yield(1) >= (1 - tolerance) % i end yields
                self.G = zeros(6,1);
self.G = [ 
(2*P(1))/self.Fy^2 + (7*M(1)^2*P(1))/(self.Fy^2*self.Mp^2);
                                                         0;
(2*M(1))/self.Mp^2 + (7*M(1)*P(1)^2)/(self.Fy^2*self.Mp^2);
                                                         0;
                                                         0;
                                                         0];
                self.ComputeKPlastic();
             elseif Yield(2) >= (1 - tolerance) % j end yields
                self.G = zeros(6,1);
                self.G  = [0;
                           0;
                           0;
                           (2*P(2))/self.Fy^2 + (7*M(2)^2*P(2))/(self.Fy^2*self.Mp^2);
                           0;
                           (2*M(2))/self.Mp^2 + (7*M(2)*P(2)^2)/(self.Fy^2*self.Mp^2)];
                self.ComputeKPlastic();
            else
                self.kp = 0;
            end
        end
        %% Calculate the plastic reduction matrix, kp
        function ComputeKPlastic(self)
            self.kp = -self.ke_local * self.G * (self.G' * self.ke_local *...
                      self.G)^-1 * self.G' * self.ke_local;
        end
    end
end