classdef CVPlanarFit
    properties
        CV
        maxCVAllowed = 2;
        biPolaeVoltageThreshold = 0.5;

    end

    methods
        function [cv cvX] = run(obj, userdata)
            % load vertices, faces and activation time from file
            %data1 = load(fname_data,'data');
            %data = data1.data;
            X = userdata.surface.triRep.X;
            faces = userdata.surface.triRep.Triangulation;
            biAct = userdata.surface.act_bip(:,1);
            biVolt = userdata.surface.act_bip(:,2);
            % tringulation is just a container for the mesh
            % vertexAttachment tneigh{i} is the list of triangles
            % attached to the i-th vertex (patch).
            tr = triangulation(faces,X);
            vn = tr.vertexNormal();
            tneigh = tr.vertexAttachments();

            vgrad = zeros(size(X,1),3);
            vnorm = zeros(size(X,1),1);
            mdelta = zeros(size(X,1),1);

            act = biAct;
            cvX = zeros(size(X,1),1);
            CVgrad = zeros(size(X,1),3);
            % iterate over all the vertices to compute the local CV
            % from the activation time. It may happen that the activation
            % is flat over the patch, so that CV is infinite. This is probably
            % noise in the data, that we smooth out by enlarging the patch to
            % its next level, until we get a reasonable approximation of CV
            CVmax = obj.maxCVAllowed;
            itermax = 1;
            vthres = obj.biPolaeVoltageThreshold;
            warning('off','all');
            for vidx = 1:size(X,1)
                % skip non-valid points
                if biVolt(vidx)<vthres || isnan(biVolt(vidx))
                    cvX(vidx) = nan;
                    continue;
                end
                % local patch and connectivity
                lver = vidx;
                lface = [];
                count = 0;
                while true
                    count = count + 1;
                    textra = [];
                    for myv = lver
                        textra = [textra; tneigh{myv}];
                    end
                    lface = unique([textra, lface]);
                    lconn = tr(lface,:);
                    lver = unique(lconn);

                    lact = act(lver)';
                    lxyz = X(lver,:);
                    % boundary vertices
                    ltri = triangulation(lconn,X);
                    bver = unique(ltri.freeBoundary());

                    % estimate CV by 1/|grad(act)| of the linear least-squares
                    % approximation of the local data
                    Xmat = [ones(length(lact),1),lxyz];
                    coef = Xmat\lact';

                    % projector on the surface: v - (v.n)n
                    PP = eye(3)-vn(vidx,:)'*vn(vidx,:);
                    cv = 1/norm(coef(2:end));
                    cvgrad = cv; %slow*cv^2;
                    fprintf('%4d %10.4f mm/ms\n',vidx,cv);

                    if count<=1, continue; end;

                    cvold = cvX(vidx);
                    cvX(vidx) = cv;
                    CVgrad(vidx,:) = cvgrad;
                    %if abs(cv-cvold) < 0.1*abs(cv) && cv < CVmax, break; end;
                    if cv < CVmax, break; end;
                    if count>=itermax
                        cvX(vidx) = nan;
                        break;
                    end

                end

                if 0
                    figure;
                    trisurf(lconn,X(:,1),X(:,2),X(:,3),act,'FaceColor','interp');
                    caxis([min(lact),max(lact)]);
                    colorbar;
                    hold on;
                    plot3(X(bver,1),X(bver,2),X(bver,3),'r*');
                    plot3(X(vidx,1),X(vidx,2),X(vidx,3),'k*');
                    break
                end

                continue

                if 1
                    % estimate based on formula min[ 1/cv*|y-x| + act(y) ]
                    DD = [];
                    DA = [];
                    PP = [];
                    p0 = X(vidx,:);
                    for ii = 1:size(bver,1)
                        vi = bver(ii);
                        vnext = bver(mod(ii,size(bver,1))+1);
                        pp = X(vi,:);
                        ss = linspace(0,1,40);
                        for s = ss
                            ps = (1-s)*pp + s*X(vnext,:);
                            as = (1-s)*act(vi) + s*act(vnext);
                            dx = norm(ps-p0);
                            dt = as-act(vidx);
                            DD = [DD, dx];
                            DA = [DA, as];
                            PP = [PP; ps];
                        end
                        % min [ 1/cv*|x-y| + u(y) ]
                    end
                    J = fminbnd(@(cc) abs(act(vidx) - min(cc*DD + DA)),-20.0,20.0);
                    disp([J,1/J]);
                    %ezplot(@(cc) abs(act(vidx)-min(cc*DD + DA)),[-20,20]);
                    %ezplot(@(cc) norm(act(vidx) + cc*DD - DA),[1,20]);
                    CC = linspace(-10,-10,1000);
                    JJ = arrayfun(@(cv) norm(act(vidx) + cv*DD - DA,2),CC);
                    [~,idx] = min(JJ);
                    disp(CC(idx));
                    disp( (PP(idx,:)-p0)/norm(PP(idx,:)-p0) );
                    figure;
                    plot(JJ);
                    return
                    plot(1:length(DD),DD);
                    hold on;
                    plot(1:length(DD),repmat(act(vidx),1,length(DD)));
                end

            end

           
            %figure

            hist(cvX,100)
            % save(fname_out,'CVxyz','bvolt','act','faces','xyz','data');
        end
         
    end
end