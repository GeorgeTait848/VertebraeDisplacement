classdef Disc
    properties
        DampingConst {mustBeNumeric}
        SpringConst {mustBeNumeric}
     
    end
    methods
        function disc = Disc(damp, spr)
            disc.DampingConst = damp;
            disc.SpringConst = spr;
        end
    end
end