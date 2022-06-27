classdef CaseLSq
    % Implements cases 1-6 as per thesis
    methods(Static)
        
        function [tc_fit, A_fit, phi0_fit, CV, residual] = Case1(times, thetacirc, radcath)
            xvec=cos(thetacirc);
            yvec=sin(thetacirc);
            timevector=times;
            
            X1 = ones(length(timevector),3);
            X1(:,2) = xvec;
            X1(:,3) = yvec;
            yt = timevector';
            beta = X1\yt;
            yhat = beta(1)+ beta(2).*X1(:,2)+beta(3).*X1(:,3);
            
            tc_fit=beta(1);
            A_fit=sqrt(beta(2)^2+beta(3)^2);
            phi0_fit=atan2(beta(3),beta(2));
            
            CV=radcath/A_fit;
            
            figure
            scatter(1:length(yt), yt, 'b')
            hold on
            scatter(1:length(yhat), yhat, 'r')
            title('Case 1 (circular catheter, planar wavefront')
            
            
            residual=sum((yt-yhat).^2);
        end
        
        function [tc_fit, A_fit,  phi0_fit, CV, residual] = Case2(times, theta, involutecircrad)
            xvec=involutecircrad.*(cos(theta)+theta.*sin(theta));
            yvec=involutecircrad.*(sin(theta)-theta.*cos(theta));
            timevector=times;
            
            X1 = ones(length(timevector),3);
            X1(:,2) = xvec;
            X1(:,3) = yvec;
            yt = timevector';
            beta = X1\yt;
            yhat = beta(1)+ beta(2).*X1(:,2)+beta(3).*X1(:,3);
            
            tc_fit=beta(1);
            AoverR_fit=sqrt(beta(2)^2+beta(3)^2);
            A_fit=AoverR_fit;
            phi0_fit=atan2(beta(3),beta(2));
            
            CV=1/AoverR_fit;
            
            figure
            scatter(1:length(yt), yt, 'b')
            hold on
            scatter(1:length(yhat), yhat, 'r')
            title('Case 2 (spiral catheter, planar wavefront')
            
            residual=sum((yt-yhat).^2);
        end
        
        function [A_fit, phi0_fit, D_fit, CV, residual] = Case3(times, thetacirc, radcath, x0)
            
            global X Y Z timev
            X=radcath^2;
            Y=2.*radcath.*cos(thetacirc);
            Z=2.*radcath.*sin(thetacirc);
            timev=times(:);
            X=X(:); Y=Y(:); Z=Z(:);
            
            beta2=lsqnonlin(@circfit,x0);   % Invoke optimizer
            beta2=real(beta2);
            A_fit=beta2(2);
            G_fit=sqrt(beta2(3)^2+beta2(4)^2);
            %D_fit=G_fit-radcath;
            D_fit=G_fit;
            phi0_fit=atan2(beta2(4),beta2(3));
            
            timessub=beta2(1) + beta2(2).*sqrt(X + (beta2(3)^2 + beta2(4)^2) +beta2(3).*Y + beta2(4).*Z);
            
            
            CV=radcath/A_fit;
            
            figure
            scatter(1:length(timev), timev, 'b')
            hold on
            scatter(1:length(timessub), timessub, 'r')
            title('Case 3 (circular catheter, circular wavefront')
            
            residual=sum((timev-timessub).^2);
        end
        
        
        function [A_fit, phi0_fit, Y_fit, CV, residual] = Case4(times, theta, involutecircrad, x0)
            
            global X Y Z timev
            X=involutecircrad^2.*(1+theta.^2);
            Y=2.*involutecircrad.*sqrt(1+theta.^2).*cos(theta-asin(involutecircrad.*theta./sqrt(involutecircrad^2.*(1+theta.^2))));
            Z=2.*involutecircrad.*sqrt(1+theta.^2).*sin(theta-asin(involutecircrad.*theta./sqrt(involutecircrad^2.*(1+theta.^2))));
            timev=times(:);
            X=X(:); Y=Y(:); Z=Z(:);
            
            beta2=lsqnonlin(@circfit,x0);   % Invoke optimizer
            beta2=real(beta2);
            %Tcombo_fit=beta2(1);
            A_fit=beta2(2);
            Y_fit=sqrt(beta2(3)^2+beta2(4)^2);
            phi0_fit=atan2(beta2(4),beta2(3));
            
            timessub=beta2(1) + beta2(2).*sqrt(X + (beta2(3)^2 + beta2(4)^2) +beta2(3).*Y + beta2(4).*Z);
            
            
            CV=1/A_fit;
            
            figure
            scatter(1:length(timev), timev, 'b')
            hold on
            scatter(1:length(timessub), timessub, 'r')
            title('Case 4 (spiral catheter, circular wavefront')
            
            residual=sum((timev-timessub).^2);
        end
        
        
        
        
        
        function [T_fit, A_fit, phi0_fit, CV, residual, beta, lowerCV] = Case5(times, xpoints, ypoints)
            [~, pos]=min(times);
            x_0=xpoints(pos);
            y_0=ypoints(pos);
            
            timevector=times;
            xvec=xpoints-x_0; %define relative to earliest activation point, check this?
            yvec=ypoints-y_0;
            
            X1 = ones(length(timevector),3);
            X1(:,2) = xvec;
            X1(:,3) = yvec;
            yt = timevector';
            beta = X1\yt;
            yhat = beta(1)+ beta(2).*X1(:,2)+beta(3).*X1(:,3);
            
            T_fit=beta(1);
            A_fit=sqrt(beta(2)^2+beta(3)^2);
            phi0_fit=atan2(beta(3),beta(2));
            
            zparameter=cos(phi0_fit).*(xvec) + sin(phi0_fit).*(yvec);
            
            [tmax, postmax]=max(yhat);
            [tmin, postmin]=min(yhat);
            
            
            
            CV=(zparameter(postmax)-zparameter(postmin))/((tmax-tmin));  %cm/ms*1000 = cm/s
            
%                         figure
%                         scatter(1:length(yt), yt, 'b')
%                         hold on
%                         scatter(1:length(yhat), yhat, 'r')
%                         title('Case 5 (arbitrary catheter, planar wavefront')
            
            residual=sum((yt-yhat).^2);
            
            [dmax, postmax]=max(zparameter);
            [dmin, postmin]=min(zparameter);
            D=(dmax-dmin);
            lowerCV=D/[max(yt)-min(yt)];
            
        end
        
        function [T_fit, A_fit, phi0_fit, D_fit, CV, beta2, residual] = Case6(times, xpoints, ypoints, x0)
            [~, pos]=min(times);
            x_0=xpoints(pos);
            y_0=ypoints(pos);
            
            global X Y Z timev
            X=(ypoints-y_0).^2 + (xpoints-x_0).^2;
            Y=2.*(xpoints-x_0);
            Z=2.*(ypoints-y_0);
            timev=times(:);
            X=X(:); Y=Y(:); Z=Z(:);
            opts = optimoptions(@lsqnonlin,'SpecifyObjectiveGradient',true, 'TolFun', 10^(-1000));
            [beta2,resnorm,residual,exitflag,output]=lsqnonlin(@circfit,x0);   % Invoke optimizer
            
            T_fit=beta2(1);
            A_fit=beta2(2);
            D_fit=sqrt(beta2(3)^2+beta2(4)^2);
            phi0_fit=atan2(beta2(4),beta2(3));%+pi;
            
            timessub=beta2(1) + beta2(2).*sqrt(X + (beta2(3)^2 + beta2(4)^2) +beta2(3).*Y + beta2(4).*Z);
            
            zparameter=sqrt(X + (beta2(3)^2 + beta2(4)^2) +beta2(3).*Y + beta2(4).*Z);
            
            [tmax, postmax]=max(timessub);
            [tmin, postmin]=min(timessub);
            
            CV=(zparameter(postmax)-zparameter(postmin))/((tmax-tmin));
            
%                         figure
%                         scatter(1:length(timev), timev, 'b')
%                         hold on
%                         scatter(1:length(timessub), timessub, 'r')
%                         title('Case 6 (arbitrary catheter, circular wavefront)')
%                         legend('actual times', 'fit times')
%                         xlabel('Electrogram number')
%                         ylabel('Relative activation time (ms)')
            
            residual=sum((timev-timessub).^2);
        end
        
        
        function [beta2, timessub] = Case8(times, xpoints, ypoints, x0)
            [~, pos]=min(times);
            x_0=xpoints(pos);
            y_0=ypoints(pos);
            
            global X Y timev
            X=xpoints;
            Y=ypoints;
            timev=times(:);
            X=X(:); Y=Y(:);
            
            beta2=lsqnonlin(@spiralfit,x0);   % Invoke optimizer
            
            
            
            timessub=beta2(1) + beta2(2).*sqrt((X-beta2(4)).^2+(Y-beta2(5)).^2)+beta2(3).*atan2((Y-beta2(5)),(X-beta2(4)));
            
            
            figure
            scatter(1:length(timev), timev, 'b')
            hold on
            scatter(1:length(timessub), timessub, 'r')
            title('Case 8 (arbitrary catheter, circular wavefront)')
            legend('actual times', 'fit times')
            xlabel('Electrogram number')
            ylabel('Relative activation time (ms)')
            
            residual=sum((timev-timessub).^2);
        end
        
        
        function [T_fit, A_fit, phi0_fit, D_fit, CV,Aniso_fit, beta2, residual] = CaseEllipse(times, xpoints, ypoints, x0)
            [~, pos]=min(times);
            x_0=xpoints(pos);
            y_0=ypoints(pos);
            
            global X1 X2 X3 X4 timev
            X1=2.*(xpoints-x_0);
            X2=2.*(ypoints-y_0);
            X3=(xpoints-x_0).^2;
            X4=(ypoints-y_0).^2;
            timev=times(:);
            X1=X1(:); X2=X2(:); X3=X3(:); X4=X4(:);
            
            beta2=lsqnonlin(@ellipsefitCV,x0);   % Invoke optimizer
            
            T_fit=beta2(1);
            A_fit=beta2(2);
            D_fit=sqrt(beta2(3)^2+beta2(4)^2);
            phi0_fit=atan2(beta2(4),beta2(3));%+pi;
            Aniso_fit=beta2(5);
            
            timessub=beta2(1) + beta2(2).*sqrt((beta2(3)^2 + beta2(4)^2) + X3 +beta2(3).*X1 + beta2(4)*beta2(5).*X2 +beta2(5)^2*X4);
            
            
            zparameter=sqrt((beta2(3)^2 + beta2(4)^2) + X3 +beta2(3).*X1 + beta2(4)*beta2(5).*X2 +beta2(5)^2*X4);
            
            [tmax, postmax]=max(timessub);
            [tmin, postmin]=min(timessub);
            
            CV=(zparameter(postmax)-zparameter(postmin))/((tmax-tmin));
            
            %             figure
            %             scatter(1:length(timev), timev, 'b')
            %             hold on
            %             scatter(1:length(timessub), timessub, 'r')
            %             title('Case 6 (arbitrary catheter, circular wavefront)')
            %             legend('actual times', 'fit times')
            %             xlabel('Electrogram number')
            %             ylabel('Relative activation time (ms)')
            
            residual=sum((timev-timessub).^2);
        end
        
        

        
        
        function [T_fit, A_fit, phi0_fit, D_fit, CV,Aniso_fit,Angle,Rot1,beta2, residual] = CaseEllipseAngleRot2(times, xpoints, ypoints, x0, b5);
            [~, pos]=min(times);
            x_0=xpoints(pos);
            y_0=ypoints(pos);
            
            global X Y timev Xmin Ymin
            X=(xpoints-x_0);
            Y=(ypoints-y_0);
            timev=times(:);
            X=X(:); Y=Y(:);
            options = optimoptions('lsqnonlin','MaxIterations', 10^10);
            if b5<1
%             lb=[-100 0.1 -100 -100 1 -1];
%             ub=[100 30 100 100 10 1];
            
%             
            lb=[-150 0.5 -100 -100 1 -1];
            lb=[-150 2/3 -50 -50 0 -1];
            ub=[150 5 50 50 3 1];
            
%             
            else
%             lb=[-100 0.1 -100 -100 0 -1];
%             ub=[100 30 100 100 1 1];        
%             
            lb=[-150 0.5 -100 -100 0 -1];
            lb=[-150 2/3 -50 -50 0 -1];
            ub=[0 5 100 100 1 1];
            ub=[150 5 50 50 3 1];
             end
            beta2=lsqnonlin(@ellipsefitCVangleRotMatNS,x0, lb,ub);   % Invoke optimizer
            T_fit=beta2(1);
            A_fit=beta2(2);
            D_fit=sqrt(beta2(3)^2+beta2(4)^2);
            phi0_fit=atan2(beta2(4),beta2(3));
            Aniso_fit=beta2(5);
            Angle=beta2(6);
            
            Rot1=asin(beta2(6));
            
            CV=1/beta2(2);
            
            timessub=beta2(1) + beta2(2).*sqrt(beta2(3)^2+beta2(4)^2 + 2.*beta2(3).*(sqrt(1-beta2(6)^2)).*X - 2.*beta2(3).*beta2(6).*Y+2.*beta2(4).*beta2(5).*beta2(6).*X+2.*beta2(4).*beta2(5).*(sqrt(1-beta2(6)^2)).*Y +((sqrt(1-beta2(6)^2).*X-beta2(6).*Y)).^2+(beta2(5).*beta2(6).*X+beta2(5).*(sqrt(1-beta2(6)^2)).*Y).^2); 
            
            residual=sum((timev-timessub).^2);
        end
      
        
            function [T_fit, A_fit, phi0_fit, D_fit, CV,Aniso_fit,Angle,Rot1,beta2, residual] = CaseEllipseAngleRot2M(times, xpoints, ypoints, x0);
            [~, pos]=min(times);
            x_0=xpoints(pos);
            y_0=ypoints(pos);
            
            global X Y timev Xmin Ymin
            X=(xpoints-x_0);
            Y=(ypoints-y_0);
            timev=times(:);
            X=X(:); Y=Y(:);
            options = optimoptions('lsqnonlin','MaxIterations', 10^10);
            lb=[-100 0.1 -100 -100 0 -1];
            ub=[100 30 100 100 10 1];
            beta2=lsqnonlin(@ellipsefitCVangleRotMatN,x0, lb,ub);   % Invoke optimizer
            T_fit=beta2(1);
            A_fit=beta2(2);
            D_fit=sqrt(beta2(3)^2+beta2(4)^2);
            phi0_fit=atan2(beta2(4),beta2(3));
            Aniso_fit=beta2(5);
            Angle=beta2(6);
            
            Rot1=asin(beta2(6));
            
            CV=1/beta2(2);
                        
            residual=0;
        end
        
    end
    
end

