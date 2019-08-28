%change to xlsread
%For example, A = sym('a',[1,3]) creates a row vector A = [a1,a2,a3].
% check for local stifness matrix
% take Iy and Iz alao from table
%% first we read the necessary data

% InM is member information matrix that has member no.(1) , frst_joint(2),
% sec_joint(3), Length(4), Area(5), depth(6), width(7), Ix(8), E(9), G(10)
M = input('Input Member file name:');
InM = xlsread(M);      %change to xlsread
% Rearranging the 2nd and 3rd collumn to have smaller joint number in
% to have smaller joint no. in 2nd collumn
for l=1:size(InM,1)
    if InM(l,3) < InM(l,2)
        temp = InM(l,3);
        InM(l,3) = InM(l,2);
        InM(l,2) = temp;
    end
end
%no. of restrained joints 
k = input('Number of Joint restrained:(2 for plane and 20 for big frame)');
% InJ is the Joint information Matrix that contain Co-ordinates, forces and momenrts of the joints.
J = input('Input Joint file name:');
InJ = xlsread(J);
Kg = zeros(6*size(InJ,1),6*size(InJ,1));

%% Making the global stiffness matrix for each member

for l=1:size(InM,1)
    %% Making the local stiffness matrix
    b = InM(l,7);
    d = InM(l,6);
    E = InM(l,11);   % Elasticity
    % Ix and Iy have interchanged formula because of taking y in horizontal
    % direction
    G=InM(l,12); J=InM(l,8); Iy=InM(1,9); Iz=InM(1,10);%% take this also from table
    L = InM(l,4);      % Length of member 
    A = InM(l,5);      % Area of cross section
    % Allocating some variables to save them in memory and prevent calculation everytime
    U = A*E/L; V = 12*Iz*E/L^3; W = 12*Iy*E/L^3;  %%     dobara dekhna h Iy aur Iz bcoz Y taken horizontal while in notebook it is vertical
    X = 6*Iz*E/L^2; Y = 6*Iy*E/L^2; Z = G*J/L;
    % Taking all elements zero at first as max. elements are zeros
    Km(:,:,l) = zeros(12,12);
    % Putting the values in the upper triangular part
    Km(1,1,l) = U; Km(1,7,l) = -U; 
    Km(2,2,l) = V; Km(2,6,l) = X; Km(2,8,l) = -V; Km(2,12,l) = X;
    Km(3,3,l) = W; Km(3,5,l) = -Y; Km(3,9,l) = -W; Km(3,11,l) = -Y;
    Km(4,4,l) = Z; Km(4,10,l) = -Z;
    Km(5,5,l) = 4*E*Iy/L; Km(5,9,l) = Y; Km(5,12,l) = 2*E*Iy/L;
    Km(6,6,l) = 4*E*Iz/L; Km(6,8,l) = -X; Km(5,12,l) = 2*E*Iz/L;
    Km(7,7,l) = U;
    Km(8,8,l) = V; Km(8,12,l) = -X;
    Km(9,9,l) = W; Km(9,11,l) = Y;
    Km(10,10,l) = Z;
    Km(11,11,l) = 4*E*Iy/L;
    Km(12,12,l) = 4*E*Iz/L;
    % Making the symetric part by transforming processes
    P = Km(:,:,l) - Km(:,:,l).*eye(12);
    Km(:,:,l) = Km(:,:,l) + P';
    
    %% Making the transformation matrix
    % Extracting near joint as nj and far joint fj
    nj = InM(l,2); fj = InM(l,3);
    % Finding the axial x-vector of local member
    VectorX = [InJ(fj,2)-InJ(nj,2) InJ(fj,3)-InJ(nj,3) InJ(fj,4)-InJ(nj,4)];
    % Converting to unit matrix
    VectorX = VectorX./norm(VectorX);
    % Saving cos of thetas as xlambdas in X, Y, Z directions
    xlambda = VectorX';
    % if x in direction of X then dt is identity 
    if xlambda(1) == 1 
        dt = eye(12);
        % if x opp. of X then dt as follows and similarly others
    elseif xlambda(1) == -1
        D = [-1 0 0;0 -1 0;0 0 1];
        dt = zeros(12,12);
        dt(1:3,1:3) = D; dt(4:6,4:6) = D; dt(7:9,7:9) = D; dt(10:12,10:12) = D;
    elseif xlambda(2) == 1
        D = [0 1 0;-1 0 0;0 0 1];
        dt = zeros(12,12);
        dt(1:3,1:3) = D; dt(4:6,4:6) = D; dt(7:9,7:9) = D; dt(10:12,10:12) = D;
    elseif xlambda(2) == -1
        D = [0 -1 0;1 0 0;0 0 1];
        dt = zeros(12,12);
        dt(1:3,1:3) = D; dt(4:6,4:6) = D; dt(7:9,7:9) = D; dt(10:12,10:12) = D;
    elseif xlambda(3) == 1
        D = [0 0 1;0 1 0;-1 0 0];
        dt = zeros(12,12);
        dt(1:3,1:3) = D; dt(4:6,4:6) = D; dt(7:9,7:9) = D; dt(10:12,10:12) = D;
    elseif xlambda(3) == -1
        D = [0 0 -1;0 1 0;1 0 0];
        dt = zeros(12,12);
        dt(1:3,1:3) = D; dt(4:6,4:6) = D; dt(7:9,7:9) = D; dt(10:12,10:12) = D;
        %Special cases for strut
    elseif xlambda(1) == 0 && xlambda(2) > 0
        Km(:,:,l) = zeros(12,12);
        Km(1,1,l) = U; Km(1,7,l) = -U;
        Km(7,1,l) = -U; Km(7,7,l) = U;
        D = [xlambda'; -1 0 0; 0 -xlambda(3) xlambda(2)];
        dt = zeros(12,12);
        dt(1:3,1:3) = D; dt(4:6,4:6) = D; dt(7:9,7:9) = D; dt(10:12,10:12) = D;
    elseif xlambda(1) == 0 && xlambda(2) < 0
        Km(:,:,l) = zeros(12,12);
        Km(1,1,l) = U; Km(1,7,l) = -U;
        Km(7,1,l) = -U; Km(7,7,l) = U;
        D = [xlambda'; 1 0 0; 0 xlambda(3) -xlambda(2)];
        dt = zeros(12,12);
        dt(1:3,1:3) = D; dt(4:6,4:6) = D; dt(7:9,7:9) = D; dt(10:12,10:12) = D;
    elseif xlambda(2) == 0 && xlambda(1) > 0
        Km(:,:,l) = zeros(12,12);
        Km(1,1,l) = U; Km(1,7,l) = -U;
        Km(7,1,l) = -U; Km(7,7,l) = U;
        D = [xlambda'; 0 1 0;-xlambda(3) 0 xlambda(1)];
        dt = zeros(12,12);
        dt(1:3,1:3) = D; dt(4:6,4:6) = D; dt(7:9,7:9) = D; dt(10:12,10:12) = D;
    elseif xlambda(2) == 0 && xlambda(1) < 0
        Km(:,:,l) = zeros(12,12);
        Km(1,1,l) = U; Km(1,7,l) = -U;
        Km(7,1,l) = -U; Km(7,7,l) = U;
        D = [xlambda'; 0 -1 0;xlambda(3) 0 -xlambda(1)];
        dt = zeros(12,12);
        dt(1:3,1:3) = D; dt(4:6,4:6) = D; dt(7:9,7:9) = D; dt(10:12,10:12) = D;
    end
    Dt(:,:,l) = dt;
    % The global matrix of each member is ready by mulltiplying with dt in
    % following manner
    K(:,:,l) = Dt(:,:,l)'*Km(:,:,l)*Dt(:,:,l);
    
    %Consider the 12x12 matrix as four big 6x6 matrix then we see that only
    %the two diagonal matrix element will add others will simply come inn
    %the matrix
    Kg(nj*6-5:nj*6,nj*6-5:nj*6) = Kg(nj*6-5:nj*6,nj*6-5:nj*6) + K(1:6,1:6,l);
    Kg(fj*6-5:fj*6,fj*6-5:fj*6) = Kg(fj*6-5:fj*6,fj*6-5:fj*6) + K(7:12,7:12,l);
    Kg(fj*6-5:fj*6,nj*6-5:nj*6) = K(7:12,1:6,l);
    Kg(nj*6-5:nj*6,fj*6-5:fj*6) = K(1:6,7:12,l);
    
    
    
end


%% Another way to create K global but upper one described is better
% for j=1:size(InJ,1)         %add condition of being equal joint make loop for all joint combination
%     for b=1:size(InJ,1)
%         if j==b
%             i = find(InM(:,2) == j);
%             for a = 1:size(i)
%                 Kg(j*6-5:j*6,j*6-5:j*6) = Kg(j*6-5:j*6,j*6-5:j*6) + K(1:6,1:6,i(a));
%             end
%             i = find(InM(:,3) == j);
%             for a = 1:size(i)
%                 Kg(j*6-5:j*6,j*6-5:j*6) = Kg(j*6-5:j*6,j*6-5:j*6) + K(7:12,7:12,i(a));
%             end
%         elseif j > b
%             i = find(InM(:,2) == b);
%             f = find(InM(:,3) == j);
%             for q = 1:size(i)
%                 for m = 1:size(f)
%                     if i(q) == f(m)                        
%                         Kg(j*6-5:j*6,b*6-5:b*6) = K(7:12,1:6,i(q));                        
%                         Kg(b*6-5:b*6,j*6-5:j*6) = K(1:6,7:12,i(q));
%                     end
%                 end
%             end
%         end
%     end
% end
% Getting all the forces 
z = InJ(1:end-k,5:10)';
%converting the forces into vectors
Pc = z(:);
%getting unknown displacements by using inverse
dis = (Kg(1:end-6*k,1:end-6*k))\Pc;
%net displacement vector
d = [dis(:); zeros(6*k,1)];
% x,y,z disp and rotations
Xd = reshape(d, 6, size(InJ,1))';
%net forces and reaction vector
P = Kg*d;
% the values of all forces and reaction in order
fr = reshape(P,6,size(InJ,1))';

for l=1:size(InM,1)
     % Extracting near joint as nj and far joint fj
    nj = InM(l,2); fj = InM(l,3);
    as = [Xd(nj,:) Xd(fj,:)]';
    f = K(:,:,l)*Dt(:,:,l)*as(:);
end