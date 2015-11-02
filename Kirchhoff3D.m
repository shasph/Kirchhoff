classdef Kirchhoff3D < handle
    
    properties
        V
        N
        F
        numPoints=0;
        numFaces=0;
                
    end
    
    methods
        function self = Kirchhoff3D(matname, t)
            
            load(matname);
            
          %  self.V = [V(:,1)*t V(:,2)/2 V(:,3)/4]; %for coin
          %  self.V = [V(:,1)*t V(:,2) V(:,3)];
            self.V = V;  
            self.F = F;
            self.N = Kirchhoff3D.renorm(self.V, self.F);            
            
            self.numPoints = length(V);
            self.numFaces = length(F);

        end
        
        function [ totalArea ] = surf_area(self)
            totalArea = Kirchhoff3D.area(self.V, self.F);
        end
        
        function [K] = compute(self, offset)
            
            C = Kirchhoff3D.face_center(self.V, self.F);

            FL = motion_flux(self);
            
            S = self.V-offset*self.N;
            
            SL = Kirchhoff3D.single_layer(S, C);
            
            %DL = Kirchhoff3D.double_layer(S, C, Kirchhoff3D.face_normal(self.V, self.F));
             
            M = Kirchhoff3D.solid_angle(S, self.V, self.F);
            
            strengths = M \ FL;
            
            phi = self.numPoints*SL*strengths;
            
            Q = one_point_quadrature(self);
            
            K = Q*phi;
            
        end

        function [MF] = motion_flux(self)
            MF = zeros(self.numFaces, 6);
            MF(:, 1:3) = Kirchhoff3D.angular_vector(self.V, self.F);
            MF(:, 4:6) = Kirchhoff3D.area_vector(self.V, self.F);
        end
        
        function [ Q ] = one_point_quadrature(self)
            n = self.numFaces;
            areas = Kirchhoff3D.triangle_area(self.V, self.F);
            VF = Kirchhoff3D.face_center(self.V, self.F);
            NF = Kirchhoff3D.face_normal(self.V, self.F);
            CR = cross(VF, NF);
            Q = zeros(6,n);
            Q(1,:) = areas.*CR(:,1);
            Q(2,:) = areas.*CR(:,2);
            Q(3,:) = areas.*CR(:,3);
            Q(4,:) = areas.*NF(:,1);
            Q(5,:) = areas.*NF(:,2);
            Q(6,:) = areas.*NF(:,3);
        end

    end
    
    methods (Static)
        
        function [ AV ] = angular_vector( pts, tris )

            n = length(tris);

            p1 = pts(tris(:, 1), :);
            p2 = pts(tris(:, 2), :);
            p3 = pts(tris(:, 3), :);

            pp1 = dot(p1', p1');
            pp2 = dot(p2', p2');
            pp3 = dot(p3', p3');

            pp12 = dot(p1', p2');
            pp23 = dot(p3', p2');
            pp31 = dot(p1', p3');

            
            p12 = (pp1+pp2+pp12)';
            p23 = (pp2+pp3+pp23)';
            p31 = (pp3+pp1+pp31)';
            
            p2p1 = p2-p1;
            p3p2 = p3-p2;
            p1p3 = p1-p3;
            
            t1(:,3) = p12.*p2p1(:,3);
            t2(:,3) = p23.*p3p2(:,3);
            t3(:,3) = p31.*p1p3(:,3);
            t1(:,1) = p12.*p2p1(:,1);
            t2(:,1) = p23.*p3p2(:,1);
            t3(:,1) = p31.*p1p3(:,1);
            t1(:,2) = p12.*p2p1(:,2);
            t2(:,2) = p23.*p3p2(:,2);
            t3(:,2) = p31.*p1p3(:,2);
            
            AV = -(t1+t2+t3)/6;
                        
        end
        
        function [ area ] = area(pts, tris)
            area = sum(Kirchhoff3D.triangle_area(pts, tris));
        end
        
        function [ areas ] = triangle_area(pts, tris)
            AV = Kirchhoff3D.area_vector( pts, tris );
            areas = sqrt(dot(AV', AV'))';
        end
        
        function [ AV ] = area_vector( pts, tris )
            p1 = pts(tris(:, 1), :);
            p2 = pts(tris(:, 2), :);
            p3 = pts(tris(:, 3), :);
            AV = 0.5*(cross(p1, p2)+cross(p2, p3)+cross(p3,p1)); 
        end
        function [ AV ] = face_normal( pts, tris )
            AV = Kirchhoff3D.area_vector( pts, tris );
            areas = Kirchhoff3D.triangle_area(pts, tris);
            AV(:,1) = AV(:,1)./areas;
            AV(:,2) = AV(:,2)./areas;
            AV(:,3) = AV(:,3)./areas;            
        end
                
        function [ SA ] = solid_angle( src, trg, tris )

            n = length(src);
            m = length(tris);

            SA = zeros(m, n);

            for i=1:n
                sp(:,1) = src(i,1)*ones(m,1);
                sp(:,2) = src(i,2)*ones(m,1);
                sp(:,3) = src(i,3)*ones(m,1);
                p1 = trg(tris(:,1),:)-sp;
                p2 = trg(tris(:,2),:)-sp;
                p3 = trg(tris(:,3),:)-sp;
                
                cr = cross(p2, p3);
                det = dot(p1', cr');

                l1 = sqrt(dot(p1',p1'));
                l2 = sqrt(dot(p2',p2'));
                l3 = sqrt(dot(p3',p3'));

                den = l1.*l2.*l3+l1.*dot(p2', p3')+l2.*dot(p3', p1')+l3.*dot(p1', p2');

                SA(:,i) = 2.0*n*atan2(det, den);
            end
        end
                
        function [ SL ] = single_layer( s, t )

            m = size(t, 1);
            n = size(s, 1);

            tx = t(:,1); ty = t(:,2); tz = t(:,3);
            sx = s(:,1); sy = s(:,2); sz = s(:,3);

            x = tx*ones(1,n) - ones(m,1)*sx';
            y = ty*ones(1,n) - ones(m,1)*sy';
            z = tz*ones(1,n) - ones(m,1)*sz';

            rr = (x.*x + y.*y + z.*z);
            
            %fprintf('l_e=%4.2f', n

            SL = 1./sqrt(rr);
        end
        
        function [ DL ] = double_layer(s, t, nt)

            m = size(t, 1);
            n = size(s, 1);

            tx = t(:,1); ty = t(:,2); tz = t(:,3);
            sx = s(:,1); sy = s(:,2); sz = s(:,3);

            x = tx*ones(1,n) - ones(m,1)*sx';
            y = ty*ones(1,n) - ones(m,1)*sy';
            z = tz*ones(1,n) - ones(m,1)*sz';

            nx = nt(:,1)*ones(1,n); ny = nt(:,2)*ones(1,n); nz = nt(:,3)*ones(1,n);

            rr = (x.*x + y.*y + z.*z);
            r3 = sqrt(rr).*rr;
            
            en = x.*nx + y.*ny + z.*nz;
            
            DL = en./r3;
        end
        
        function [ CP ] = face_center(v, f)
            CP = (v(f(:,1),:)+v(f(:,2),:)+v(f(:,3),:))/3;
        end
        
%        function [ NT ] = face_normal(v, f)
%             NT = cross(V(F(:,2),:)-V(F(:,1),:), V(F(:,3),:)-V(F(:,1),:));
%             NT = normr(fn);
%         end
        
        %calculate the vercter normal based on
        %http://www.flipcode.com/archives/Vertex_Normals.shtml
        function [N] = renorm(v,f)
            nv = length(v);
            nf = length(f);
            
            vn = zeros(nv,3);
            
            fn = cross(v(f(:,2),:)-v(f(:,1),:), v(f(:,3),:)-v(f(:,1),:));
            fn = normr(fn);
            
            for k=1:nv
                for p=1:nf
                    if (k==f(p,1))||(k==f(p,2))||k==(f(p,3))
                        vn(k,:)= vn(k,:)+fn(p,:);
                    end
                end
                vn(k,:)=vn(k,:)/norm(vn(k,:));
            end
            N = vn;
        end
        
    end
end

