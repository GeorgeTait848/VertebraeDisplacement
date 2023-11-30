classdef Spine
    properties
        BottomVertebrae {mustBeVertebrae}
        Array {mustBeVector}
        TopVertebrae {mustBeVertebrae}
    end
    methods
        function spine = Spine(BottomVertebrae, Array, TopVertebrae)
            spine.BottomVertebrae = BottomVertebrae;
            spine.Array = Array;
            spine.TopVertebrae = TopVertebrae;
        end
    end
end