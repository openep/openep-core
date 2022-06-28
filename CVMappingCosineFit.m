classdef CVMappingCosineFit
    %function [cv, cvX, interpCv] = doCvMapping_CosineFit( userdata, int )
    % DOCVMAPPING_COSINEFIT Calculates conduction velocities using
    % triangulation of electrode activation times
    %
    % Usage:
    %   [cv, cvX, interpCv] = doCvMapping_CosineFit( userdata, int )
    % Where:
    %   userdata - see importcarto_mem
    %   int      - see openEpDataInterpolator.m
    %   cv       - the calculated conduction velocity data, in m/s
    %   cvX      - the Cartesian co-ordinates at which conduction velocity data
    %              has been calculated. size(cvX) = [length(cv), 3].
    %   interpCv - conduction velocity data interpolated across the surface of
    %              the shell.
    %              size(interpCv) = [length(userdata.surface.triRep.X), 1].
    %
    % DOCVMAPPING_COSINEFIT Calculatess conduction velocities using
    % cosine fit technique
    %
    % Author:
    % SPDX-License-Identifier: Apache-2.0
    %
    % Modifications -
    %
    % Info on Code Testing:
    % -----------------------------------s----------------------------
    %
    % ---------------------------------------------------------------
    %
    % ---------------------------------------------------------------
    % code
    % ---------------------------------------------------------------

    % first perform global interpolation using cosine fit to
    % calculate conduction velocities

    % TODO

    % accept only those conduction velocity values in proximity to electrodes?
    % TODO. We should in some way limit cv and cvX only to values that are
    % likely to be real; i.e. in close proximity; or at; mapping points.

    % now do interpolation, using the interpolator specified, so we have a full
    % dataset at each of the mesh nodes.
    %vtx = getVertices(userdata, 'used', false);

    properties
        plot = true;
    end

    methods

        function [cv, cv_triangle] = run(obj, userdata)
           
        newTri=userdata.surface.triRep.Triangulation;
        newX=userdata.surface.triRep.X; %mm

        LAT_Bi=userdata.surface.act_bip;
        LAT=LAT_Bi(:,1); %local activation times
        BipolarVoltage=LAT_Bi(:,2);

        figure
        p=patch('Vertices', newX, 'Faces', newTri);
        set(p,'facecolor','blue','edgecolor','black');
        daspect([1 1 1]);
        view([180 0]); axis tight; grid on;
        camlight; lighting gouraud;
        alpha(.25)
        hold on
        scatter3(newX(:,1), newX(:,2), newX(:,3), 10)

        UsedVerts=unique(newTri(:));
        AllVerts=1:length(newX);
        Unused=setdiff(AllVerts, UsedVerts);

        NewVerts=newX;
        NewVerts(Unused, :)=[];

        newTriRelabel=newTri; newTriRelabel(:)=nan;

        for kk=1:length(newTriRelabel)
            F1=newTri(kk,1); F1_L=find(F1==UsedVerts);
            F2=newTri(kk,2); F2_L=find(F2==UsedVerts);
            F3=newTri(kk,3); F3_L=find(F3==UsedVerts);
            newTriRelabel(kk,1)=F1_L;
            newTriRelabel(kk,2)=F2_L;
            newTriRelabel(kk,3)=F3_L;
        end

        figure
        p=patch('Vertices', NewVerts, 'Faces',newTriRelabel);
        set(p,'facecolor','blue','edgecolor','black');
        daspect([1 1 1]);
        view([180 0]); axis tight; grid on;
        camlight; lighting gouraud;
        alpha(.25)

        newX=NewVerts;
        newTri=newTriRelabel;
        LAT(Unused)=[];

        %plot LATs interpolated/extrapolated by the system

        figure
        p=patch('Vertices', newX, 'Faces', newTri);
        set(p,'facecolor','blue','edgecolor','black');
        daspect([1 1 1]);
        view([180 0]); axis tight; grid on;
        camlight; lighting gouraud;
        alpha(.25)
        hold on
        scatter3(newX(:,1), newX(:, 2), newX(:, 3), 50, LAT, 'filled')
        colorbar
        colormap('jet')
        
        %alternatively work with recorded points only
        Locations=userdata.electric.egmSurfX;
        Timings=userdata.electric.annotations.mapAnnot;
        Tags=userdata.electric.tags;
        ToDelAbl=[];
        for inds=1:length(Tags)
            TT=Tags{inds};
            if length(TT)>0
                ToDelAbl=[ToDelAbl inds];

            end
        end

        Exclude=find(Timings==-8000);
        Exclude=unique([Exclude; ToDelAbl']);

        Locations(Exclude,:)=[];
        Timings(Exclude)=[];
        Timings=Timings-min(Timings);

        close all
        figure
        p=patch('Vertices', newX, 'Faces', newTri);
        set(p,'facecolor','blue','edgecolor','black');
        daspect([1 1 1]);
        view([180 0]); axis tight; grid on;
        camlight; lighting gouraud;
        alpha(.25)
        hold on
        scatter3(Locations(:,1), Locations(:, 2), Locations(:, 3), 50, Timings, 'filled')
        colorbar
        colormap('jet')

        connects=newTri;
        verts=newX;

        clearvars -except connects verts Locations Timings LoopOver

        Septum_XYZ=Locations;
        Septum_LAT=Timings;

        FMPX=(verts(connects(:,1), 1)+verts(connects(:,2), 1)+verts(connects(:,3), 1))./3;
        FMPY=(verts(connects(:,1), 2)+verts(connects(:,2), 2)+verts(connects(:,3), 2))./3;
        FMPZ=(verts(connects(:,1), 3)+verts(connects(:,2), 3)+verts(connects(:,3), 3))./3;


        CP=zeros(1, length(Septum_LAT));
        for kk=1:length(Septum_LAT)
            D=(Septum_XYZ(kk,1)-FMPX).^2+(Septum_XYZ(kk,2)-FMPY).^2+(Septum_XYZ(kk,3)-FMPZ).^2;
            [~, loc]=min(D);
            CP(kk)=loc;
        end

        CPLAT=CP;

        surfaceX(:,1)=FMPX(CP);
        surfaceX(:,2)=FMPY(CP);
        surfaceX(:,3)=FMPZ(CP);

        figure
        p=patch('Vertices',verts, 'Faces',connects);
        set(p,'facecolor','blue','edgecolor','none');
        daspect([1 1 1]);
        view(3); axis tight; grid on;
        camlight; lighting gouraud;
        alpha(.5)
        hold on
        scatter3(surfaceX(:,1), surfaceX(:,2), surfaceX(:,3), 50,Septum_LAT, 'filled')
        colormap('jet')
        colorbar


        cd( '/home/ali/Documents/MATLAB/numerical-tours-master/matlab')
        getd = @(p)path(p,path);

        getd('toolbox_signal/');
        getd('toolbox_general/');
        getd('toolbox_graph/');
        getd('toolbox_wavelet_meshes/');


        heartTriangles=connects;
        heartVertices=verts;

        LAT=Septum_LAT;

        figure
        p=patch('Vertices', heartVertices, 'Faces', heartTriangles);
        set(p,'facecolor','blue','edgecolor','black');
        daspect([1 1 1]);
        view([180 0]); axis tight; grid on;
        camlight; lighting gouraud;
        alpha(.25)
        hold on
        scatter3(surfaceX(:,1), surfaceX(:, 2), surfaceX(:, 3), 50, LAT, 'filled')
        colorbar
        colormap('jet')

        rovingX=Septum_XYZ;

        figure
        p=patch('Vertices', heartVertices, 'Faces', heartTriangles);
        set(p,'facecolor','blue','edgecolor','none');
        daspect([1 1 1]);
        view([180 0]); axis tight; grid on;
        camlight; lighting gouraud;
        alpha(.25)
        hold on
        scatter3(rovingX(:,1), rovingX(:, 2), rovingX(:, 3), 50, LAT, 'filled')
        colorbar
        colormap('jet')

        %select points in 2cm sphere, geodesic flattening

        Radius=6; % 8 changed
        AllCVStore=zeros(length(surfaceX),9);

        List=1:length(surfaceX);

        ActStore=LAT;


        FMP=[FMPX FMPY FMPZ];

        RovP=zeros(length(rovingX), 1);

        for kk=1:length(rovingX)
            D=(rovingX(kk,1)-FMP(:,1)).^2+(rovingX(kk,2)-FMP(:,2)).^2+(rovingX(kk,3)-FMP(:,3)).^2;
            [~, loc]=min(D);
            RovP(kk)=loc;
        end

        ToFix=RovP;
        length(ToFix)

        figure
        p=patch('Vertices', heartVertices, 'Faces', heartTriangles);
        set(p,'facecolor','blue','edgecolor','none');
        daspect([1 1 1]);
        view([180 0]); axis tight; grid on;
        camlight; lighting gouraud;
        alpha(.25)
        hold on
        scatter3(FMP(RovP,1), FMP(RovP, 2), FMP(RovP, 3), 50, LAT, 'filled')
        colorbar
        colormap('jet')

        close all

        for index=1:length(ToFix)
            cd('/home/ali/Desktop/Codes')

            %for each point to check, find 10 closest points within a time window of the closest
            %point. Add in a modulus function for working across a cycle length

            clearvars -except LoopOver rovingX surfaceX ToFix FMP ActStore Radius heartVertices heartTriangles AllCVStore VertInd List VI CPLAT index LAT
           

            index

            VertInd=(ToFix(index));

            CentrePoint=(ToFix(index));

            Dmat=(heartVertices(:,1)-FMP(CentrePoint,1)).^2+(heartVertices(:,2)-FMP(CentrePoint,2)).^2+(heartVertices(:,3)-FMP(CentrePoint,3)).^2;
            Dmat=sqrt(Dmat);
            Points=find(Dmat<Radius);

           


            %
            vertices=heartVertices;
            faces=heartTriangles;

            %take all vertices in a 1cm x 1cm box around it (2 x 2)
            Point=heartVertices(Points, :);
            EGRAMX=Point(:, 1);
            EGRAMY=Point(:, 2);
            EGRAMZ=Point(:, 3);
            delta=0.1; delta=0.05; delta=delta*10;
            minX=min(EGRAMX)-delta;
            maxX=max(EGRAMX)+delta;
            minY=min(EGRAMY)-delta;
            maxY=max(EGRAMY)+delta;
            minZ=min(EGRAMZ)-delta;
            maxZ=max(EGRAMZ)+delta;
            FindVertX=(minX<vertices(:,1)).*(vertices(:,1)<maxX);
            FindVertY=(minY<vertices(:,2)).*(vertices(:,2)<maxY);
            FindVertZ=(minZ<vertices(:,3)).*(vertices(:,3)<maxZ);
            FindVert=FindVertX.*FindVertY.*FindVertZ;
            FindVertList=find(FindVert);

            
            %subselection
            FaceList=[];
            for kk=1:length(FindVertList)
                F1=find(faces(:,1)==FindVertList(kk));
                F2=find(faces(:,2)==FindVertList(kk));
                F3=find(faces(:,3)==FindVertList(kk));
                FL=[F1; F2; F3];
                FaceList=[FaceList; FL];
                FaceList=unique(FaceList);
            end


            NewFaces=faces(FaceList, :);
            NewList=unique(NewFaces(:));
            NewVerts=vertices(NewList, :);

            FFound=find(FaceList==VertInd);

          

            %find times in this window
           

            FoundT=intersect(CPLAT, FaceList);
            LocsT=zeros(1, length(FoundT));
            for ijk=1:length(LocsT)
                FF=find(FoundT(ijk)==CPLAT);
                LocsT(ijk)=min(FF);
            end





           
            for kk=1:length(NewFaces)
                NewFaces(kk,1)=find(NewList==NewFaces(kk,1));
                NewFaces(kk,2)=find(NewList==NewFaces(kk,2));
                NewFaces(kk,3)=find(NewList==NewFaces(kk,3));
            end

            MP3DX=(NewVerts(NewFaces(:,1), 1)+NewVerts(NewFaces(:,2), 1)+NewVerts(NewFaces(:,3), 1))./3;
            MP3DY=(NewVerts(NewFaces(:,1), 2)+NewVerts(NewFaces(:,2), 2)+NewVerts(NewFaces(:,3), 2))./3;
            MP3DZ=(NewVerts(NewFaces(:,1), 3)+NewVerts(NewFaces(:,2), 3)+NewVerts(NewFaces(:,3), 3))./3;

            

            vertex=NewVerts'; Ffaces=NewFaces';
            options.niter=10000;
            W=ones(length(vertex), 1); I=8;
            [U,err,Usvg] = perform_geodesic_iterative(vertex, Ffaces, W, I, options);

          
            [path,vlist,plist] = compute_geodesic_mesh(U, vertex, Ffaces, 1, options);

            
            %
            Dstore=cell(1, length(vertex));
        
            kk=1;
            while (max(U)<100)&&(kk<(length(vertex)+1))
                I=kk;
                [U,err,Usvg] = perform_geodesic_iterative(vertex, Ffaces, W, I, options);
                Dstore{kk}=U;
                kk=kk+1;
            end

            if max(U)<100

                Dmat=zeros(length(Dstore), length(Dstore));
                for kk = 1:length(Dstore)
                    Dmat(kk, 1:length(Dstore))=(Dstore{kk});
                end
                geod=Dmat;

                %make a symmetric matrix
                Dup=triu(Dmat);
                DupF=Dup+Dup';

                Dmin=DupF;


                if length(find(isinf(Dmin)))>0
                else
                    [Y,stress] = mdscale(Dmin,2, 'criterion','sammon');

                    newx=Y(:, 1);
                    newy=Y(:, 2);
                 
                    MP2DX=(newx(NewFaces(:,1))+newx(NewFaces(:,2))+newx(NewFaces(:,3)))./3;
                    MP2DY=(newy(NewFaces(:,1))+newy(NewFaces(:,2))+newy(NewFaces(:,3)))./3;

                    %ActTime=ActStore(NewList(SelectP)); any associated activation times?

                    %relate new mid-points to old faces
                    NewListOfFaces=FaceList;

                    ActTimeList=zeros(1, length(FaceList));

                    for kk=1:length(FaceList)
                        Found=min(find(FaceList(kk)==CPLAT));
                        if length(Found)>0
                            ActTimeList(kk)=Found;
                        end
                    end



                    ListToUse=ActTimeList(ActTimeList>0);

                    if length(ListToUse)<3
                    else
                        xpo=MP2DX(find(ActTimeList>0));
                        ypo=MP2DY(find(ActTimeList>0));

                        ActivationTimes=LAT(ListToUse);

                       
                        %
                        Locations3D=[MP3DX(find(ActTimeList>0)) MP3DY(find(ActTimeList>0)) MP3DZ(find(ActTimeList>0))];
                        
                        %baseline fit:
                        rowtouse=xpo; coltouse=ypo;
                        ActTime=ActivationTimes;


                        xpointsmm=coltouse.*1;
                        ypointsmm=rowtouse.*1;
                        times=ActTime';
                        times=times-min(times);


                        [T_fit5, A_fit5, phi0_fit5, CV5, residual5, beta] = CaseLSq.Case5(times, xpointsmm, ypointsmm);
                        [val, pos]=min(times);
                        x_0=xpointsmm(pos);
                        y_0=ypointsmm(pos);
                        lengtharrow=2;
                        u1=lengtharrow*cos(phi0_fit5);
                        v1=lengtharrow*sin(phi0_fit5);
                        CVplanar=CV5;

                        beta10=0;
                        beta20=A_fit5;
                        Dguess=1;
                        beta30=Dguess*cos(phi0_fit5);
                        beta40=Dguess*sin(phi0_fit5);
                        x0 = [beta10 beta20 beta30 beta40];

                        xvec=xpointsmm-x_0; %define relative to earliest activation point, check this?
                        yvec=ypointsmm-y_0;

                        yhat = beta(1)+ beta(2).*xvec+beta(3).*yvec;
                        yt=times;
                        MeanRes=sum((yt(:)-yhat(:)).^2)/length(yt);
                        ResidualsPW=(yt(:)-yhat(:)).^2;
                        MedianRes=mean(ResidualsPW);

                        %remove points with residual > 2x mean residual.
                        ToRemoveRes=find(ResidualsPW>2*MedianRes);



                        times(ToRemoveRes)=[];
                        xpointsmm(ToRemoveRes)=[];
                        ypointsmm(ToRemoveRes)=[];
                        ActivationTimes(ToRemoveRes)=[];

                        Locations3D(ToRemoveRes, :)=[];

                       

                        times=times-min(times);

                        %now take a smaller subset

                        %central point?


                        [T_fit5, A_fit5, phi0_fit5, CV5, residual5, beta] = CaseLSq.Case5(times, xpointsmm, ypointsmm);
                        [val, pos]=min(times);
                        x_0=xpointsmm(pos);
                        y_0=ypointsmm(pos);
                        lengtharrow=2;
                        u1=lengtharrow*cos(phi0_fit5);
                        v1=lengtharrow*sin(phi0_fit5);
                        CVplanar=CV5;

                        beta10=0;
                        beta20=A_fit5;
                        Dguess=1;
                        beta30=Dguess*cos(phi0_fit5);
                        beta40=Dguess*sin(phi0_fit5);
                        x0 = [beta10 beta20 beta30 beta40];

                        xvec=xpointsmm-x_0; %define relative to earliest activation point, check this?
                        yvec=ypointsmm-y_0;

                        yhat = beta(1)+ beta(2).*xvec+beta(3).*yvec;
                        yt=times;
                        MeanRes=sum((yt(:)-yhat(:)).^2)/length(yt);


                        [T_fit, A_fit, phi0_fit, D_fit, CV, beta, residual]=CaseLSq.Case6(times, xpointsmm, ypointsmm, x0);

                        CVcircular=CV;
                        [val, pos]=min(times);
                        x_0=xpointsmm(pos);
                        y_0=ypointsmm(pos);
                        x1=x_0; y1=y_0;
                        lengtharrow=D_fit;
                        x2=x1-lengtharrow*cos(phi0_fit);
                        y2=y1-lengtharrow*sin(phi0_fit);



                        MeanRes=residual5/length(times);





                        if length(xpointsmm)>4



                            [T_fit5, A_fit5, phi0_fit5, CV5, residual5, beta] = CaseLSq.Case5(times, xpointsmm, ypointsmm);
                            [val, pos]=min(times);
                            x_0=xpointsmm(pos);
                            y_0=ypointsmm(pos);
                            lengtharrow=2;
                            u1=lengtharrow*cos(phi0_fit5);
                            v1=lengtharrow*sin(phi0_fit5);
                            CVplanar=CV5;

                            beta10=0;
                            beta20=A_fit5;
                            Dguess=1;
                            beta30=Dguess*cos(phi0_fit5);
                            beta40=Dguess*sin(phi0_fit5);
                            x0 = [beta10 beta20 beta30 beta40];


                            [T_fit, A_fit, phi0_fit, D_fit, CV, beta, residual]=CaseLSq.Case6(times, xpointsmm, ypointsmm, x0);

                            CVcircular=CV;
                            [val, pos]=min(times);
                            x_0=xpointsmm(pos);
                            y_0=ypointsmm(pos);
                            x1=x_0; y1=y_0;
                            lengtharrow=D_fit;
                            x2=x1-lengtharrow*cos(phi0_fit);
                            y2=y1-lengtharrow*sin(phi0_fit);

                         

                            CX=mean(xpointsmm); CY=mean(ypointsmm);
                            DOrigin=(x_0-CX)^2+(y_0-CY)^2;
                            Dsource=(x2-CX)^2+(y2-CY)^2;

                            SourceInCatheter=Dsource<DOrigin;

                            %add in source location

                            S1=x_0-D_fit*cos(phi0_fit);
                            S2=y_0-D_fit*sin(phi0_fit);


                            VV=[newy newx zeros(size(newx))];

                           
                            %express in barycentric coords

                            % %find if a point lies in a triangle,
                            % http://www.blackpawn.com/texts/pointinpoly/

                            P=[MP2DY(FFound).*1+u1./1 MP2DX(FFound).*1+v1./1]; %0.5?
                            vert2d=VV(:, 1:2).*1;
                            F2=NewFaces;

                            v0 = vert2d(F2(FFound,3), :) - vert2d(F2(FFound,1), :);
                            v1 = vert2d(F2(FFound,2), :) - vert2d(F2(FFound,1), :);
                            v2 = P - vert2d(F2(FFound,1), :);


                            % Compute dot products
                            dot00 = dot(v0, v0);
                            dot01 = dot(v0, v1);
                            dot02 = dot(v0, v2);
                            dot11 = dot(v1, v1);
                            dot12 = dot(v1, v2);

                            % Compute barycentric coordinates
                            invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
                            u = (dot11 * dot02 - dot01 * dot12) * invDenom;
                            v = (dot00 * dot12 - dot01 * dot02) * invDenom;

                            % Check if point is in triangle
                            cond=1-((u > 0) && (v > 0) && (u + v < 1));


                            Location=vert2d(F2(FFound, 3),:).*u+vert2d(F2(FFound, 2),:).*v+vert2d(F2(FFound, 1),:).*(1-u-v);
                            

                            Location3d=NewVerts(F2(FFound, 3),:).*1.*u+NewVerts(F2(FFound, 2),:).*1.*v+NewVerts(F2(FFound, 1),:).*1.*(1-u-v);
                           
                            u=Location3d(1)-MP3DX(FFound).*1;
                            v=Location3d(2)-MP3DY(FFound).*1;
                            w=Location3d(3)-MP3DZ(FFound).*1;
                            
                           

                            A1=find(ActTimeList>0);
                            
                            u=Location3d(1)-MP3DX(FFound).*1;
                            v=Location3d(2)-MP3DY(FFound).*1;
                            w=Location3d(3)-MP3DZ(FFound).*1;
                            

                            %                 %modify st choose points based on geodesic distances

                            starttime=1;
                            CVstore(starttime)=CVplanar;
                            resstore(starttime)=residual5;
                            Ustore(starttime)=u;
                            Vstore(starttime)=v;
                            Wstore(starttime)=w;
                            NoPointsC(starttime)=length(times);
                           

                            CVstoreC(starttime)=CVcircular;
                            DstoreC(starttime)=D_fit;
                            resstoreC(starttime)=residual;




                            CVstore(CVstore==0)=[]; resstore(resstore==0)=[]; CVstore=CVstore'; resstore=resstore';
                            Ustore(Ustore==0)=[]; Vstore(Vstore==0)=[]; Ustore=Ustore'; Vstore=Vstore';
                            Wstore(Wstore==0)=[]; Wstore=Wstore';
                            NoPointsC(NoPointsC==0)=[];

                            %add in CV circular, d0 and residual circular

                            CVstoreC(CVstoreC==0)=[]; resstoreC(resstoreC==0)=[]; CVstoreC=CVstoreC'; resstoreC=resstoreC';
                            DstoreC(DstoreC==0)=[]; DstoreC=DstoreC';


                            CVstoreNoPoints=CVstore;
                            ResstoreNoPoints=resstore;
                            NostoreNoPoints=NoPointsC;
                            UstoreNoPoints=Ustore;
                            VstoreNoPoints=Vstore;
                            WstoreNoPoints=Wstore;
                            CVstoreNoPointsC=CVstoreC;
                            ResstoreNoPointsC=resstoreC;
                            DstoreNoPoints=DstoreC;

                        else
                            CVstoreNoPoints=nan;
                            ResstoreNoPoints=nan;
                            NostoreNoPoints=nan;
                            UstoreNoPoints=nan;
                            VstoreNoPoints=nan;
                            WstoreNoPoints=nan;

                            CVstoreNoPointsC=nan;
                            ResstoreNoPointsC=nan;
                            DstoreNoPoints=nan;
                        end



                        clearvars CVstore resstore NoPointsC Ustore Vstore Wstore CVstoreC resstoreC DstoreC


                        CollateCV=(CVstoreNoPoints(:));
                        CollateRes=ResstoreNoPoints(:);
                        CollateCV(CollateRes==0)=[];
                        CollateRes(CollateRes==0)=[];
                        CollateU=UstoreNoPoints(:);
                        CollateU(CollateU==0)=[];
                        CollateV=VstoreNoPoints(:);
                        CollateV(CollateV==0)=[];
                        CollateW=WstoreNoPoints(:);
                        CollateW(CollateW==0)=[];

                        CollateCVC=CVstoreNoPointsC(:);
                        CollateCVC(CollateCVC==0)=[];

                        CollateResC=ResstoreNoPointsC(:);
                        CollateResC(CollateResC==0)=[];

                        CollateD=DstoreNoPoints(:);
                        CollateD(CollateD==0)=[];

                        CollateNo=NostoreNoPoints(:);
                        CollateNo(CollateNo==0)=[];


                        if length(CollateCV)==length(CollateRes)
                            CollateAll=[CollateCV CollateRes CollateU CollateV CollateW CollateCVC CollateResC CollateD CollateNo];


                            ToKeep=find((CollateAll(:,2)<200000));
                            CollateAll2=CollateAll;
                            CollateAll2=CollateAll2(ToKeep,:);
                            CollateAll2(CollateAll2(:,1)>5.0, :)=[]; %1.5
                            %CollateAll2=unique(CollateAll2);

                            AllCVStore(index, 1:length(CollateAll2))=CollateAll2;
                        end

                    end

                end
               
            else
                AllCVStore(index, :)=nan;
            end



        end


        save(strcat('CV_data', num2str(LoopOver), '.mat'), '-v7.3')
        end
    end
end



    %end