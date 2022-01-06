function [interp_values, d_interp_values]=rbf_interpolator(values,x_values,x_all,rbfoptions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uses a modified version of the "Scattered data interpolation and approximation using RBFs"
% (https://uk.mathworks.com/matlabcentral/fileexchange/10056-scattered-data-interpolation-and-approximation-using-radial-base-functions)
% functions on the MATLAB Fileexchange to interpolate LATs from measurement
% points to triangulated surface, and gradients of LATs in order to compute
% velocities:
% 
% Based on the scheme in: -------------------------------------------------
% M. Masé, F. Ravelli, Automatic reconstruction of activation and velocity 
% maps from electro-anatomic data by radial basis functions,
% in: Annual International Conference of the IEEE Engineering in Medicine and Biology Society (EMBC),
% 2010, pp. 2608?2611, doi:10.1109/IEMBS.2010.5626616.
% -------------------------------------------------------------------------
%
% Inputs:
% in1: scalar field (1,n) of samples where n is number of samples
% in2: locations of sample points (d,n) where d is dimension
% in3: locations of points to interpolate to (d,N)
% in4: Radial basis function (RBF) type ('multiquadric' works well)
% in5: scaling constant for RBF (testing shows 1 works well, but it's problem dependent)
%     
% Outputs:
% out1: interpolated scalar field (1,N)
% out2: gradient of interpolated scalar field (d,N)
% [Outputs 3 to 5 are specific to cardiac mapping]
% out3: wave velocity (d,N)
% out4: wave direction (d,N) (the unit-vector field)
% out5: velocity magnitude, or speed (1,N)
%
% Notes: this function indicates the quality of the interpolation, if this
% much greater than 1e-9 then play with the RBF scale parameter ('const')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Modified from 'RBFConductionVeloctiy.m' by Christopher O'Shea 2021
% If lats = x_lats, then function will just fit rbf field to the points

% Some error-checking...could probably be improved:
rbf=rbfoptions.basisFunction;
const=rbfoptions.shapeParameter;
optimisation=rbfoptions.doOptimisation;


fun=@(c)rbfcheck2(x_values',values',rbf,c)
if optimisation == true
const = fmincon(fun,const);
end
op = rbfcreate(x_values', values', 'RBFFunction',rbf,'RBFConstant',const);
assignin('base','op',op)
[interp_values, d_interp_values] = rbfinterp(x_all', op);
end

function [err]=rbfcheck2(x_lats,lats,rbf,const)
options=rbfcreate(x_lats, lats, 'RBFFunction',rbf,'RBFConstant',const);  
s = rbfinterp(x_lats, options);
assignin('base','s',s)
%err=max(abs(s-lats));
err=sqrt(mean((s - lats).^2));
end