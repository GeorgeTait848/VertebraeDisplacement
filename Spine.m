classdef Spine
    properties
        Vertebrae (1,:) Vertebrate
        Discs (1,:) Disc
        len (1,1) {mustBeInteger, mustBePositive} = 1
        
    end
    methods

        function spine = Spine(Vertebrae, Discs)

            if length(Discs) ~= length(Vertebrae) - 1
                error("Must have exactly one fewer disk than vertebrae")
            end

            spine.Vertebrae = Vertebrae;
            spine.Discs = Discs;
            spine.len = length(Vertebrae);
        end
        

        function y = middleVertAndAdjacentStates(self, i)

            arguments
                self Spine
                i (1,1) {mustBeInteger, mustBePositive}
            end
            
            y = zeros(3,2);

            for j=1:3
                
                idxInSpine = i + j - 2;
                y(j,:) = self.Vertebrae(idxInSpine).getState();

            end

        end

        function y = bottomVertAndAdjacentStates(self)
            
            y = zeros(2,2);

            for j=1:2
                y(j,:) = self.Vertebrae(j).getState();
            end

        end

        function y = topVertAndAdjacentStates(self)
            
            y = zeros(2,2);
            
            for i=1:2
                idxInSpine = self.len - i + 1;
                y(i,:) = self.Vertebrae(idxInSpine).getState();
            end

        end

        function dydt = middleVertDerivs(self, i)
            
            arguments
                self Spine
                i (1,1) {mustBeInteger, mustBePositive}
        
            end
            
            currDisc = self.Discs(i);
            k = currDisc.SpringConst;
            gamma = currDisc.DampingConst;
            y = self.middleVertAndAdjacentStates(i);

            dydt = zeros(1,2);
            dydt(1) = y(2,2);
            dydt(2) = -k*(2*y(2,1) - y(1,1) - y(3,1)) - gamma*(2*y(2,2) - y(1,2) - y(3,2));
            

        end

        function dydt = topVertDerivs(self)

            currDisc = self.Discs(self.len - 1);
            k = currDisc.SpringConst;
            gamma = currDisc.DampingConst;
            y = self.topVertAndAdjacentStates;

            dydt = zeros(1,2);
            dydt(1) = y(2,2);
            dydt(2) = -k*(y(2,1) - y(1,1)) - gamma*(y(2,2) - y(1,2));
           
        
        end

        function dydt = bottomVertDerivs(self, t)

            arguments
                self Spine
                t (1,1) {mustBePositive}
            end
            currDisc = self.Discs(1);
            k = currDisc.SpringConst;
            gamma = currDisc.DampingConst;
            y = self.topVertAndAdjacentStates;

            %TODO: Find a more realistic implementation of the driving
            %force from engine vibration
            h = @(x) (5*sin((pi/2)*(x-3))) .* ((x>3) & (x<5));

            dydt = zeros(1,2);
            dydt(1) = y(2,2);
            dydt(2) = -k*(y(1,1) - y(2,1)) - gamma*(y(1,2) - y(2,2)) + h(t);
        end 


        function dydt = getFullDerivs(self, t)
            %%% Returns the velocity, acceleration of each vertebrate at
            %%% any time t, as a nx2 vector.
            %%% n being the number of vertebrae in the spine.
            %%% ith row is the velocity, acceleration of ith vertebrate
            %%% the first element in each row is velocity
            %%% second is the acceleration
            arguments
                self Spine
                t (1,1) {mustBePositive}
            end

            dydt = zeros(self.len, 2);
            dydt(1) = self.bottomVertDerivs(t);


            for i = 2:self.len-1
                dydt(i) = self.middleVertDerivs(i);
            end

            dydt(self.len) = self.topVertDerivs();

        end
    end

end