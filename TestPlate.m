function TestPlate()
    %size of plate
    a = 50;
    b = 30;
    c = 2;%thickness
    %estimated realmass according to papers
    AA = RealAdd(a,b,c);
    EE = RealEll(a,b,c);
%     output = fopen('result.txt','w');
    for n=2%1:2
        filename = sprintf('plate%d.mat', n);
        K = Kirchhoff3D(filename);
        %Scale vertices to match area of plate
        K.V = sqrt(RealArea(a,b,c)/Kirchhoff3D.area(K.V, K.F))*K.V;
        for k=3%1:5
%             offset = c*0.01*k;
            offset = 0.3;
            tic;
            KT = K.compute(offset);
            KK = diag(KT);        
            dt = toc;
            save('tensor.mat','KT');
%             fprintf('%s offset=%f [t=%f] |V|=%d |F|=%d\n', filename, offset, dt, length(K.V), length(K.F));
%             fprintf(output, '%s offset=%f [t=%f] |V|=%d |F|=%d\n', filename, offset, dt, length(K.V), length(K.F));
%             fprintf(output, 'K11=%0.12f K22=%0.12f K33=%0.12f K44=%0.12f K55=%0.12f K66=%0.12f\n', KK(1), KK(2), KK(3), KK(4), KK(5), KK(6));
        end
    end
%     fclose(output);
end
        
function [ qual_err, rel_err ] = error( K , AA)
    dd = diag(K);
    rel_err = norm(AA - dd)/norm(AA);
    ndd = norm(dd);
    dd = dd./ndd;
    nAA = norm(AA);
    AA = AA./nAA;
    da = abs(AA - dd);
    qual_err = norm(da);
end

%calculate the area of the plate
function [S] = RealArea(a,b,c)
  S = 2*(a*b + b*c + a*c);
end

%real value of added-mass
function [R] = RealAdd(a,b,c)
m11 = pi*c*c/4;
m22 = pi*b*b/4;
m33 = pi*a*a/4;
I11 = pi*(b^2-c^2)^2*a/128;
I22 = pi*(c^2-a^2)^2*b/128;
I33 = pi*(a^2-b^2)^2*c/128;

R = [I33;I22;I11;m33;m22;m11];
end

%real value of ellipsoid(similar size)
function [E] = RealEll(a,b,c)
sum = a*b*c;

A = @(x)1./(((a.^2+x).^1.5).*((b.^2+x).^0.5).*((c.^2+x).^0.5));
B = @(x)1./(((b.^2+x).^1.5).*((a.^2+x).^0.5).*((c.^2+x).^0.5));
C = @(x)1./(((c.^2+x).^1.5).*((b.^2+x).^0.5).*((a.^2+x).^0.5));

a0 = sum*quadgk(A,0,Inf); 
b0 = sum*quadgk(B,0,Inf); 
c0 = sum*quadgk(C,0,Inf); 
i0 = 4*pi*sum/15;

m11 = 4*a0*pi*sum/(3*(2-a0));
m22 = 4*b0*pi*sum/(3*(2-b0));
m33 = 4*c0*pi*sum/(3*(2-c0));
m44 = i0*inertia(b,c,b0,c0);
m55 = i0*inertia(c,a,c0,a0);
m66 = i0*inertia(a,b,a0,b0);

E = [m66; m55; m44; m33; m22; m11];
end

%for RealEll function
function I = inertia(x,y,x0,y0)
t1 = (x^2-y^2)^2*(y0-x0);
t2 = 2*(x^2-y^2) + (x^2+y^2)*(x0-y0);
if t2 == 0 
    I=0;
else
    I = t1/t2;
end
end

    
    
    




