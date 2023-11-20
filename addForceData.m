function userdata = addForceData(userdata)

time = zeros(size(userdata.rf.originaldata.ablparams.time));
force = zeros(size(userdata.rf.originaldata.ablparams.time));
axialangle = zeros(size(userdata.rf.originaldata.ablparams.time));
lateralangle = zeros(size(userdata.rf.originaldata.ablparams.time));
position = zeros(size(userdata.rf.originaldata.ablparams.time));

forcedata.time = time;
forcedata.force = force;
forcedata.axialangle = axialangle;
forcedata.lateralangle = lateralangle;
forcedata.position = position;

userdata.rf.originaldata.force = forcedata;

end