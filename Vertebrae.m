classdef Vertebrae
    properties
        Mass {mustBePositive}
        Displacement {mustBeNumeric}
    end
    methods
        function vert = Vertebrae(mass, displ)
            disc.Mass = mass;
            disc.Displacement = displ;
        end
    end
end