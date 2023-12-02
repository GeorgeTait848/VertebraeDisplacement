classdef Vertebrate
    properties

        Mass (1,1) {mustBePositive} = 1
        % Need to set default mass value to 1 as MATLAB stupidly assigns
        % default value for scalar numeric values to 0, which doesnt meet
        % the mustBePositive Requirement. 
        Displacement (1,1) {mustBeNumeric} = 0
        Velocity (1,1) {mustBeNumeric} = 0
    end
    methods
        function vert = Vertebrate(mass, displ, vel)
            vert.Mass = mass;
            vert.Displacement = displ;
            vert.Velocity = vel;
        end

        function y = getState(self)

            y = zeros(1,2);
            y(1) = self.Displacement;

        end

        function self = updateState(self, y)

            arguments 
                self Vertebrate
                y (1,2) {mustBeNumeric}
            end

            self.Displacement = y(1);
            self.Velocity = y(2);
        end
    end
end