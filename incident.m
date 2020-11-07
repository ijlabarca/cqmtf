classdef incident < handle

    properties
    end

    methods
        
        % This is an abstract base class for the plane wave and
        % point source type incident fields. The methods
        % just do the minimum... the idea is that they
        % get replaces by specific versions for the other
        % incident fields.

        %-----------------------------------------------
        % constructor
        %-----------------------------------------------

        function cof = get_coefficients(self,centre,nmax)

            cof=zeros(2*nmax+1,1);

        end


        %-----------------------------------------------
        % evaluate method
        %-----------------------------------------------

        function val = evaluate(self,points,mask)

            val=zeros(size(points));

        end

        %-----------------------------------------------
        % evaluate method
        %-----------------------------------------------

        function [dx,dy] = evaluateGradient(self,points,mask)

            error('Not implemented yet')

            dx = zeros(size(points));
            dy = zeros(size(points));

        end

    end

end
