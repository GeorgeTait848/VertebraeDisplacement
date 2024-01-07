clear; 
timeOfSim = 50;
dt = 0.002;
numVertebrae = 23;
vertMass = 1;
springConst = 1;
dampConst = 1;


vertebrae = repelem(Vertebrate(vertMass, 0, 0), numVertebrae);
discs = repelem(Disc(dampConst, springConst), numVertebrae-1);
spine = Spine(vertebrae, discs);

y = spine.sim(timeOfSim, dt);
size(y)