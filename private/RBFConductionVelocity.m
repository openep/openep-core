function [interpLATs,d_interpLATs,u,n,speed]=RBFConductionVelocity(lats,x_lats,x_all,varargin)
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

% Some error-checking...could probably be improved:
if nargin==5
    rbf=varargin{1};
    const=varargin{2};
elseif nargin==4
    rbf=varargin{1};
    const=1;
elseif nargin==3
    rbf='multiquadric';
    const=1;
else
    error('Need inputs 1-3...see help: ''help RBFConductionVelocity''');
end

op = rbfcreate(x_lats, lats, 'RBFFunction',rbf,'RBFConstant',const);
rbfcheck(op); %if this returns anything >10^-9 then modify 'const'
[interpLATs, d_interpLATs] = rbfinterp(x_all, op);
mag_df=sqrt(d_interpLATs(1,:).^2+d_interpLATs(2,:).^2+d_interpLATs(3,:).^2);
speed=1./mag_df;
n=[speed.*d_interpLATs(1,:);speed.*d_interpLATs(2,:);speed.*d_interpLATs(3,:)];
u=[speed.*n(1,:);speed.*n(2,:);speed.*n(3,:)];
end