classdef Disc
    properties
        DampingConst {mustBeNumeric}
        SpringConst {mustBeNumeric}
        Position {mustBeInteger}
    end
    methods
        function disc = Disc(damp, spr, pos)
            disc.DampingConst = damp;
            disc.SpringConst = spr;
            disc.Position = pos;
        end
    end
end