classdef solver < handle

    properties
        kwave
        incidentField
        numIncidentField
    end

    methods

        %===============================================================
        % methods that will (probably) not be overridden in child class
        %===============================================================

        %-----------------------------------------
        % constructor
        %-----------------------------------------

        function self = solver(kwave_real, kwave_imag, incidentField)

            self.kwave = kwave_real +1i*kwave_imag;

            if ~iscell(incidentField)

                self.incidentField{1} = incidentField;

            else

                self.incidentField = incidentField;

            end

            self.numIncidentField = length(self.incidentField);

        end


        function setIncidentField(self,incidentField)

            % make sure that the incident field is stored as a cell array
            % even if there is only one incident field
            if ~iscell(incidentField)

                self.incidentField{1} = incidentField;

            else

                self.incidentField = incidentField;

            end

            % set number of incident fields
            self.numIncidentField = length(self.incidentField);

        end


    end % end methods

    methods(Abstract=true)

        %===============================================================
        % these methods must be overridden in the child class
        %===============================================================

        %-----------------------------------------
        % setup
        %-----------------------------------------

        setup(self)

        %-----------------------------------------
        % solve
        %-----------------------------------------

        solve(self)


    end % end methods

end
