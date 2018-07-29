% Cesar Y. Marco & David Tse (CMDT)
% CEE 282: Programming Project
% 2d2in Node Class
% 03/18/16

classdef CMDT_Node_2d2in < RC_Node_2d1el
% This is the child class of RC_Node_2d1el.m.
% All the protected & public data properties of parent class will be 
% inherited to this class 

% Node class for 2nd order analysis of a 2-dimensional structure
    
    % Protected properties go here
    properties (Access = protected)
    end
    
    % Public methods go here
    methods (Access = public)
        % Constructor function inherits 2d1el parameters
        function self = CMDT_Node_2d2in(node_number, node_coord)
            self = self@RC_Node_2d1el(node_number, node_coord);
        end
        % Updates the x and y coordinates only
        function UpdateCoordinates(self,new_node_coord)
            self.node_coord = new_node_coord;
        end 
    end
end
