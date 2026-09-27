%Code to test ideas about fmap
%Uday, started 01 May 2020
%v2: works for upto four objects (JOSA A paper setting)
%v3: included mpme
%v4: included incident field formulation, new inc field expression
clearvars -except cfg
close all

run_obj = false; %generate objects in the sim
run_fwd = false; %run the foward solver
run_grd = false; %generate synthetic data over grid
gen_invp_spec = false; %generate inverse problem geometry
run_pts = true; %generate the rand locations and corresp fields
run_pts_mpme = false; %generate mpme points
run_inv = false; %run the inverse solver on above locations
run_ttf = false; %calculate true tangential fields
run_prd = false; %generate the prediction matrix


nobj = 4; %set the number of objects here
snr = 25; %dB
lambda = 0.2; %wavelength
mtd_id = 0; %0 - cssom; 1 - nm6
Ni = 5; %number of inc field coefficients
inc_id = 1; %0-point source;1-square source
sr = 0.3; %sampling rate

scheme = 'mpme'; %sampling scheme used for the inverse problem: 'mpme' or 'random'
if exist('cfg','var') %settings passed in by run_all.m override the defaults above
    fn = fieldnames(cfg);
    for ii = 1:numel(fn), eval([fn{ii} ' = cfg.(fn{ii});']); end
end

fname_o = strcat('geom_',num2str(nobj),'.mat'); %file name for objects
fname_f = strcat('fwd_',num2str(nobj),'.mat'); %file name for forward system
fname_p = strcat('pts_',num2str(sr),'_',num2str(nobj),'.mat'); %file name for locations
fname_g = strcat('sgrid_lam_10',num2str(nobj),'.mat'); %file name for grid sol
fname_i = strcat('invp_',num2str(nobj),'.mat');%file name for inv prob specs
fname_s = strcat('syt_',num2str(sr),'_',num2str(nobj),'.mat'); %file name for synth fields
fname_a = strcat('state_',num2str(nobj),'.mat');%file name for state equation
fname_t = strcat('trtang_',num2str(nobj),'.mat');%file name for true tang fields
fname_e = strcat('prdmtx_',num2str(nobj),'.mat');%file name for prediction matrix
fname_mpme  = strcat('pts_mpme_',num2str(sr),'_',num2str(nobj),'.mat');%file name for prediction matrix
fname_mp_s  = strcat('syt_mpme_',num2str(sr),'_',num2str(nobj),'.mat');%file name for prediction matrix

show_vis = false; vis = [];
if show_vis
    figure; vis = gca; hold; grid;
end
plot_2d = true; %plots 2D fields



%some conventions: outermost region is R0, air region is R1, first object
%is R2, and so on. The interface between Region i and i+1 is indexed by i,
%so S0 is between R0 and R1, etc. Same holds true for the normals. All
%normal vectors point into R1 (i.e. except S0, all are outward normals).

%% this creates all objects and datastructures of geometry
if run_obj
    tic
    dl_fwd = lambda/40;
    %Matrial properties
    eps = zeros(2+nobj,1);
    eps(1) = 3.7 - 2.1i; %wall
    eps(2) = 1; %air
    if nobj == 1
        eps(3) = 3.7 - 2.1i; %object 1
    elseif nobj == 2
        eps(3) = 3.7 - 2.1i; %object 1
        eps(4) = 1.7 - 1.1i; %object 2
    elseif nobj == 3
        eps(3) = 3.7 - 2.1i; %object 1
        eps(4) = 1.7 - 1.1i; %object 2
        eps(5) = 2.7 - 3.7i; %object 3
    elseif nobj == 4
        eps(3) = 3.7 - 2.1i; %object 1
        eps(4) = 1.7 - 1.1i; %object 2
        eps(5) = 2.7 - 3.7i; %object 3
        eps(6) = 1.2 - 0.1i; %object 4
    else
        error('nobj exceed maximum limit (4).');
    end
    geom.eps = eps;
    %Surfaces (list of segment coordinates and normal vectors)
    geom.surf = struct([]); %for surfaces
    geom.corners = struct([]); %for corners of bounding boxes
    geom.shape = []; %for shape of object
    %for rect shape, it stores bottom left and top right corners in BL and TR
    %for circ shape, it stores center in BL and radius in TR(1). TR(2) = 0
    %Create your objects here
    for i = 1:nobj+1
        if i==1 %Wall -- S0
            shp = 'rect'; %shape of contour
            BL = [0,0]; TR = [10,10];
            nfac = -1; %for this surface alone we want inward normal
        elseif i==2 %Object 1 -- S1
            shp = 'rect'; %shape of contour
            BL = [2.5,6]; TR = [3.5,7];
            nfac = 1;
        elseif i==3 %Object 2 -- S2
            shp = 'circ'; %shape of contour
            BL = [7,7]; TR = [0.75,0]; %center at (2,2), radius of 0.75
            nfac = 1;
        elseif i==4 %Object 3 -- S3
            shp = 'rect'; %shape of contour
            BL = [5,2.125]; TR = [7,2.875];
            nfac = 1;
        elseif i==5 %Object 4 -- S4
            shp = 'circ'; %shape of contour
            BL = [2.5,4]; TR = [1,0];
            nfac = 1;
        end
        BL = BL * lambda; TR = TR * lambda;
        [p,m,q,nv] = gen_surfgrid(BL,TR,dl_fwd,vis,shp);
        nv = nfac * nv;
        geom.surf = [geom.surf struct('p',p,'m',m,'q',q,'n',nv)];
        geom.corners = [geom.corners struct('bl',BL,'tr',TR)];
        geom.shape = [geom.shape; shp];
    end
    geom.lambda = lambda;
    geom.dl = dl_fwd;
    if inc_id == 0
        geom.source = struct('bl',[lambda/2+5*lambda -3*lambda/4+5*lambda]...
            ,'tr',[lambda/2+5*lambda -3*lambda/4+5*lambda]);
    elseif inc_id == 1
        geom.source = struct('bl',[-lambda/8+5*lambda,-lambda/8+5*lambda]...
            ,'tr',[lambda/8+5*lambda,lambda/8+5*lambda]);
    else
        error('Invalid inc_id.');
    end
    save(fname_o,'geom','lambda','dl_fwd');
    disp(['Generated ' num2str(nobj) ' object(s) and stored to disk ' fname_o]);
    toc
else
    load(fname_o);
end

%% This solves the forward problem
if run_fwd
    tic
    nreg = numel(geom.eps);%number of regions
    %count number of variables, n(i) has the nvars of surf i-1
    nvar = zeros(nreg-1,1);
    for i = 1:nreg-1
        nvar(i)= 2 * length(geom.surf(i).m(:,1));
    end
    N = sum(nvar);
    nsur = length(nvar); %no of surfaces
    %Build system matrix: e.g with 1 obj, sys is block 4x2
    %Region: testing pts : source points
    %R0: testing S0: [S0-S0] 0       = 0
    %R1: testing S0: [S0-S0] [S0-S1] = inc(S0)
    %R1: testing S1: [S1-S0] [S1-S1] = inc(S1)
    %R2: testing S1: 0       [S1-S1] = 0
    %
    %In general, cols = nsur, nrows = (nreg-1)+nsur. First term is the
    %number of regions other than R1 which contribute 1 row, and for R1 we
    %have the number of surfaces giving the same number of rows. The rows
    %from 2 to 1+nsur will be fully dense, but not 1 and 2+nsur to nrows
    %     k = 2*pi*sqrt(geom.eps(i))/lambda; %wavevector
    %order of variables is [dphi/dn_0,phi_0,dphi/dn_1,phi_1]
    %Store the "blueprint" structure of the above sys matrix
    %Format: kvec, region, testing pts, source pts, strt position in matrix
    nr = 0; nc = 0;
    bpt = struct('kve',2*pi*sqrt(geom.eps(1))/lambda,'reg',0,...
        'tst',geom.surf(1),'src',geom.surf(1),'pos',[nr nc]); %First row
    %Get ready for Rows of R1
    nr = size(geom.surf(1).m,1);
    for i=1:nsur
        nc = 0;
        for j=1:nsur
            bpt = [bpt, struct('kve',2*pi*sqrt(geom.eps(2))/lambda,'reg',1,...
                'tst',geom.surf(i),'src',geom.surf(j),'pos',[nr,nc])]; %Rows of R1
            nc = nc + 2 * size(geom.surf(j).m,1);
        end
        nr = nr + size(geom.surf(i).m,1);
    end
    %Get ready for Rows for Rx, x>1
    nc = 2 * size(geom.surf(1).m,1);
    for i=1:nsur-1 %Rows for Rx, x>1
        bpt = [bpt, struct('kve',2*pi*sqrt(geom.eps(2+i))/lambda,'reg',i+1,...
            'tst',geom.surf(1+i),'src',geom.surf(1+i),'pos',[nr,nc])];
        nr = nr + size(geom.surf(i+1).m,1);
        nc = nc + 2 * size(geom.surf(i+1).m,1);
    end
    
    %Now make the system matrix
    A = zeros(N,N);
    b = zeros(N,1);
    [w,z] = glquadrule(2); %Use an n-pt quadrature rule
    for s = 1:length(bpt)
        n_tst = size(bpt(s).tst.p,1);
        for ii=1:n_tst
            n_src = size(bpt(s).src.p,1);
            %disp(['ii ' num2str(bpt(s).pos(1)+ii)]);
            for jj=1:n_src
                %intg_green(r,p,q,n,k,reg,w,z)
                [g,dg] = intg_green(bpt(s).tst.m(ii,:),bpt(s).src.p(jj,:),...
                    bpt(s).src.q(jj,:),bpt(s).src.n(jj,:),bpt(s).kve,bpt(s).reg,w,z);
                A(bpt(s).pos(1)+ii,bpt(s).pos(2)+jj) = g;
                A(bpt(s).pos(1)+ii,bpt(s).pos(2)+jj+n_src) = -dg;
                %disp(['jj ' num2str(bpt(s).pos(2)+jj) ',' num2str(bpt(s).pos(2)+jj+n_src)]);
            end
        end
    end
    %Now the RHS vector
    rhs_ctr = size(geom.surf(1).m,1);
    for i=1:nsur
        for j=1:size(geom.surf(i).m,1)
            %assuming sources are only present in region 1
            b(rhs_ctr+j) = incfield(geom.eps(2),lambda,geom.surf(i).m(j,:),inc_id);
        end
        rhs_ctr = rhs_ctr + j;
    end
    %Solve
    x = A\b;
    save(fname_f,'A','b','x','lambda','bpt');
    disp(['Solved system eqs at lambda/' num2str(lambda/dl_fwd) ...
        ' and stored to ' fname_f]);
    toc
else
    load(fname_f);
end

%% this generates the true soln over a grid
if run_grd
    tic;
    %this generates the field on a grid
    dl_g = lambda/10;
    extents = geom.corners(1).tr - geom.corners(1).bl;
    nx = floor(extents(1)/dl_g);
    ny = floor(extents(2)/dl_g);
    x_g = dl_g/2:dl_g:(nx-dl_g/2)*dl_g;
    y_g = fliplr(dl_g/2:dl_g:(ny-dl_g/2)*dl_g);
    [X_g,Y_g] = meshgrid(x_g,y_g);
    locs_g = [];
    ind_ex = [];
    for i = 1:numel(X_g)
        r_g = [X_g(i),Y_g(i)];
        if in_air(r_g,geom,0)
            locs_g = [locs_g; r_g];
        end
    end
    A_grid = gen_matrices(geom,locs_g,false,inc_id);
    sol_grid = A_grid * x;
    sol_gridi = zeros(size(sol_grid));
    for i = 1:size(locs_g,1)
        sol_gridi(i) = incfield(geom.eps(2),geom.lambda,locs_g(i,:),inc_id);
    end
    save(fname_g,'sol_grid','sol_gridi','locs_g','dl_g')
    disp(['Saved true soln on grid to ' fname_g]);
    if plot_2d == true
        Etr = zeros(size(X_g));
        for i = 1:size(X_g,1)
            for j = 1:size(X_g,2)
                ind1 = find(locs_g(:,1)==X_g(i,j) & locs_g(:,2)==Y_g(i,j));
                if isempty(ind1)
                    Etr(i,j) = 0.01;
                else
                    Etr(i,j) = sol_grid(ind1);
                end
            end
        end
        figure;
        imagesc(abs(Etr));
        toc;
    end
    
else
    load(fname_g);
end

%% Generate or load inverse problem geometry
if gen_invp_spec
    dl_inv = lambda/5; %res to solve inv prob
    inex = 0.1 * lambda; % specify nonzero amount for inexact boundaries
    invp.shape = [];
    for i=1:length(geom.surf)
        sgn = -1;
        if i == 1
            sgn = 1; %inexact outer most surface must go inward for i=1
        end
        invp.corners(i) = geom.corners(i); %same corners for surfaces
        if geom.shape(i,:) == 'rect'
            shp_i = 'rect'; %shape in inverse problem
            invp.corners(i).bl = invp.corners(i).bl + sgn* [inex,inex];
            invp.corners(i).tr = invp.corners(i).tr - sgn* [inex,inex];
        else
            shp_i = 'rect'; %shape in inverse problem
            rad = geom.corners(i).tr(1,1);
            invp.corners(i).bl = geom.corners(i).bl + sgn*[rad+inex,rad+inex];
            invp.corners(i).tr = geom.corners(i).bl - sgn*[rad+inex,rad+inex];
        end
        %gen surfaces at dl_inv
        BL = invp.corners(i).bl;
        TR = invp.corners(i).tr;
        [p,m,q,n] = gen_surfgrid(BL,TR,dl_inv,vis,shp_i);
        if i==1
            nfac = -1; %this is to make the normal on S0 inward pointing
        else
            nfac = +1;
        end
        n = nfac*n;
        invp.surf(i) = struct('p',p,'m',m,'q',q,'n',n);
        invp.shape = [invp.shape; shp_i];
    end
    invp.dl = dl_inv;
    invp.inex = inex;
    invp.eps = geom.eps(2);
    invp.lambda = lambda;
    invp.source = geom.source;
    save(fname_i,'invp');
else
    load(fname_i);
end
is_inex = false;
if abs(invp.inex)/lambda > 1e-3
    is_inex = true;
end
if is_inex
    disp(['Working with inexact boundaries (' num2str(invp.inex/invp.lambda) ...
        '*lambda) for the inverse problem at grid lambda/' num2str(invp.lambda/invp.dl)]);
else
    disp(['Working with exact boundaries for the inverse problem at grid lambda/' ...
        num2str(invp.lambda/invp.dl)]);
end

N_o = zeros(numel(invp.surf),1);
for i = 1:length(N_o)
    N_o(i) = length(invp.surf(i).m);
end


%% generate a set of random points and store exact field for them
if run_pts
    tic;
    nmeas = round(sr*2*sum(N_o));
    excl = 0.05*lambda; %exclude this much region around each boundary
    excl_s = 0.05*lambda; %exclude this much region around existing pt
    x_L = invp.corners(1).bl(1,1) + excl;
    x_R = invp.corners(1).tr(1,1) - excl;
    y_B = invp.corners(1).bl(1,2) + excl;
    y_T = invp.corners(1).tr(1,2) - excl;
    i = 0; %counter for no of points added
    locs = zeros(nmeas,2); %store the locations here, [x,y] per row
    while i < nmeas
        rx = x_L + (x_R-x_L)*rand();
        ry = y_B + (y_T-y_B)*rand();
        %check if legit
        legit = in_air([rx,ry],invp,excl);
        %now check self distances
        for j=1:i
            if norm([rx,ry] - locs(j,:)) < excl_s
                legit = false;
            end
        end
        %Add if all is satisfied
        if legit == true
            i = i+1;
            locs(i,:) = [rx,ry];
        end
    end
    disp(['Stored ' num2str(size(locs,1)) ' locs to ' fname_p ', excl=' ...
        num2str(excl) ',' num2str(excl_s)]);
    if show_vis
        scatter(vis,locs(:,1),locs(:,2),'k+');
    end
    save(fname_p,'locs','excl','excl_s');
    
    %now generate the synthetic data for these points (total field, as for MPME)
    A_est = gen_matrices(geom,locs,false,inc_id);
    b_ff = A_est * x; %estimated noiseless scattered field at locs
    i_ff = zeros(size(b_ff));
    for i = 1:size(locs,1)
        i_ff(i) = incfield(invp.eps,invp.lambda,locs(i,:),inc_id);
    end
    b_nf = sampled_noise(b_ff+i_ff,snr);
    sig = sqrt(norm(b_nf)^2/(10^(snr/10)));
    A_est = gen_matrices(invp,locs,false,inc_id);
    save(fname_s,'A_est','b_ff','b_nf','sig','snr');
    disp(['Stored synthetic & noisy scattered fields to ' fname_s]);
    toc
elseif run_pts_mpme
    dl_g = lambda/10;
    excl = 0.1*lambda;
    extents = invp.corners(1).tr-excl - invp.corners(1).bl+excl;
    nx = floor(extents(1)/dl_g);
    ny = floor(extents(2)/dl_g);
    x_g = dl_g/2:dl_g:(nx-dl_g/2)*dl_g;
    y_g = dl_g/2:dl_g:(ny-dl_g/2)*dl_g;
    [X_g,Y_g] = meshgrid(x_g,y_g);
    locs_pr = [];
    
    for i = 1:numel(X_g)
        r_g = [X_g(i),Y_g(i)];
        if in_air(r_g,invp,excl)
            locs_pr = [locs_pr; r_g];
        end
    end
    
    Bs = gen_matrices(invp,locs_pr,false,inc_id);
    Bi = gen_incmtx(invp,locs_pr,Ni,false,inc_id);
    B = [Bs,Bi];    
    vnorm = vecnorm(B);
    B = B*diag(1./vnorm);
    [U,S,V]  =svd(B);

    n=1; xs = 0; convg = false;
    while ~convg
        sigma(n) = S(n,n);
        Res(n)  = sum(sigma.^2)/sum(diag(S).^2);
        
        if Res < 0.999
            n = n+1;
        else
            convg = true;
        end
    end
    Pr=conj(V(:,1:n));
    B1 = B*Pr;
    selected = [];
    available = 1:size(B1,1); % N = size(B1,1)
    
    %% Select the first n-1 sensing locations
    phi = [];
    P = eye(n);
    for k = 1:n-1
        ip = zeros(1,length(available));
        for i = 1:length(available)
            ip(i) = norm(P*B1(available(i),:)');
        end
        [~,maximum] = max(ip);
        phi = [phi;B1(available(maximum),:)];
        selected = [selected, available(maximum)];
        available(maximum)=[];
        R = orth(phi');
        P = eye(n) - R*R';
    end
    
    %% Determine the remaining sensing locations
    lam = -10;
    tol = 0.12;
    while (lam<tol)
        
        [V,D] = eig(phi'*phi);
        % Determine multiplicity of minimum eigenvalue
        U = V(:,1);
        for i = 2:n
            if(D(i)~=D(1))
                break;
            end
            U = [U,V(i,:)];
        end
        P = U*U';
        ip = zeros(1,length(available));
        for i = 1:length(available)
            ip(i) = norm(P*B1(available(i),:)');
        end
        [~,maximum] = max(ip);
        phi = [phi;B1(available(maximum),:)];
        selected = [selected, available(maximum)];
        if length(selected) == 500
            break;
        end
        available(maximum)=[];
        lam = D(1);
    end
    
    %% Display the selected points and the number
    disp(length(selected));
    % Grid of points at lambda/40
    nmeas = round(sr*2*sum(N_o));
    locs_mpme_all(:,1) = locs_pr(selected,1); locs_mpme_all(:,2) = locs_pr(selected,2);
    locs_mpme(:,1) = locs_mpme_all(1:nmeas,1); locs_mpme(:,2) = locs_mpme_all(1:nmeas,2);
    if show_vis
        scatter(vis,locs_mpme(:,1),locs_mpme(:,2),'k+'); 
    end
    disp(['Stored MPME data to ' fname_mpme]);
    locs = locs_mpme;
    save(fname_mpme,'locs','B','locs_pr');
    %now generate the synthetic data for these points
    A_est = gen_matrices(geom,locs,false,inc_id);
    b_ff = A_est * x; %estimated noiseless scattered field at locs
    i_ff = zeros(size(b_ff));
    for i = 1:nmeas
        i_ff(i) = incfield(invp.eps,invp.lambda,locs(i,:),inc_id);
    end
    b_nf = sampled_noise(b_ff+i_ff,snr);
    sig = sqrt(norm(b_nf)^2/(10^(snr/10)));
    A_est = gen_matrices(invp,locs,false,inc_id);
    save(fname_mp_s,'A_est','b_ff','b_nf','sig','snr');
    disp(['Stored synthetic & noisy scattered fields to ' fname_mp_s]);    
else
%     load(fname_p);
%     load(fname_s);
    load(fname_mpme);
    load(fname_mp_s);
end

%% Run the inverse problem on the set of random points
if run_inv
    tic;
    % check state error of true soln on inexact boundary
    if is_inex
        % get the true soln on inexact boundaries
        if run_ttf
            ptsoni = []; %these are the midpoints of the segments
            nsoni_1 = []; nsoni_2 = []; %pts for the normals to compute dphi/dn
            h = invp.dl/50; %factor computing one sided finite differences
            for i = 1:numel(invp.surf)
                ptsoni  = [ptsoni; invp.surf(i).m]; %take mid pt of segment
                nsoni_1 = [nsoni_1; invp.surf(i).m + h*invp.surf(i).n];
                nsoni_2 = [nsoni_2; invp.surf(i).m - h*invp.surf(i).n];
            end
            plot_inex = false;
            if plot_inex
                figure; scatter(ptsoni(:,1),ptsoni(:,2)); grid; hold;
                scatter(nsoni_1(:,1),nsoni_1(:,2),'xr')
                scatter(nsoni_2(:,1),nsoni_2(:,2),'+k'); axis equal
            end
            [A_ptsoni,~]  = gen_matrices(geom,ptsoni,false,inc_id); %must be called with geom
            [A_nsoni_1,~] = gen_matrices(geom,nsoni_1,false,inc_id); %must be called with geom
            [A_nsoni_2,~] = gen_matrices(geom,nsoni_2,false,inc_id); %must be called with geom
            for i = 1 : length(nsoni_1)
                inc_field(i,1)  = incfield(geom.eps(2),lambda,ptsoni(i,:),inc_id);
                inc_field1(i,1) = incfield(geom.eps(2),lambda,nsoni_1(i,:),inc_id);
                inc_field2(i,1) = incfield(geom.eps(2),lambda,nsoni_2(i,:),inc_id);
            end
            x_oni = A_ptsoni * x  ;
            x_noni_1 = A_nsoni_1 * x ;
            x_noni_2 = A_nsoni_2 * x ;
            for i=1 :length(x_oni) %add the incident field to get total field
                x_oni(i) = x_oni(i) + incfield(invp.eps,lambda,ptsoni(i,:),inc_id);
                x_noni_1(i) = x_noni_1(i) + incfield(invp.eps,lambda,nsoni_1(i,:),inc_id);
                x_noni_2(i) = x_noni_2(i) + incfield(invp.eps,lambda,nsoni_2(i,:),inc_id);
            end
            %now x_oni contains the phi components of the field, next dphi/dn
            xd_oni = (x_noni_1 - x_noni_2)/(2*h);
            x_full_oni = zeros(2*length(x_oni),1);
            off1 = 0; off2 = 0;
            for i=1:numel(invp.surf)
                Li = length(invp.surf(i).m);
                x_full_oni(off1+1:off1+Li) = xd_oni(off2+1:off2 +Li);
                x_full_oni(off1+Li+1:off1+2*Li) = x_oni(off2+1:off2 +Li);
                off1 = off1 + 2*Li;
                off2 = off2 + Li;
            end
            save(fname_t,'x_full_oni','ptsoni');
        else
            load(fname_t);
        end
    else
        for i = 1:numel(invp.surf)
            ptsoni  = [ptsoni; invp.surf(i).m]; %take mid pt of segment
        end
        x_full_oni = zeros(length(ptsoni),1); % change later to downsample
    end
    [A_states,b_state] = gen_matrices(invp,ptsoni,true,inc_id);
    A_statei = gen_incmtx(invp,ptsoni,Ni,true,inc_id);
    A_state = [A_states, A_statei];
    srt = norm(A_states*x_full_oni - b_state)/norm(b_state);
    disp(['State residual with true soln on inexact is ' num2str(srt)]);
    
    %build datastructures for solving inv prob
    [A_ests,~] = gen_matrices(invp,locs,false,inc_id); %A_est eqn
    A_esti = gen_incmtx(invp,locs,Ni,false,inc_id);
    A_est = [A_ests, A_esti];

%     b_inv = A_est * x_full_oni;
    toc;
    
    %solve for the tangential fields
    D = zeros(2*sum(N_o)+Ni);
    N_pr = 0;
    for i = 1:length(N_o)
        D(N_pr+1:N_pr+N_o(i),N_pr+1:N_pr+N_o(i)) = dctmtx(N_o(i));
        D(N_pr+N_o(i)+1:N_pr+2*N_o(i),N_pr+N_o(i)+1:N_pr+2*N_o(i)) = dctmtx(N_o(i));
        N_pr = N_pr + 2*N_o(i);
    end
    D(N_pr+1:end,N_pr+1:end) = dctmtx(Ni);
    
    vnorm = vecnorm([A_est;A_state]);
%     vnorm = ones(1,length(vnorm));
    A_est = A_est*diag(1./vnorm); A_state = A_state*diag(1./vnorm);
    save(fname_a,'A_state','b_state','A_est','D','vnorm');

    if mtd_id==0
        disp('Solving with CSSOM');
        [U,S,V] = svd(A_est);
        %now implement morozov's principle
        fac = 1;
        p=1; xs = 0; convg = false;
        while ~convg
            xs = xs + V(:,p) * (U(:,p)'*b_nf)/S(p,p);
            res = norm(A_est*xs-b_nf);
            if res > fac*sig
                p = p+1;
            else
                convg = true;
            end
        end
        %now implement SOM and L1 minimization
        Vs = V(:,1:p);
        fac1 = 0.6; fac_s = 0.1;
%         ze1 = zeros(2*sum(N_o),1); ze2 = zeros(Ni,1);
        I1 = eye(2*sum(N_o)+Ni); I1(2*sum(N_o)+1:end,2*sum(N_o)+1:end) = zeros(Ni);
        I2 = eye(2*sum(N_o)+Ni); I2(1:2*sum(N_o),1:2*sum(N_o)) = zeros(2*sum(N_o)); 
        Q = eye(size(A_est,2))-Vs*Vs'; %the projector matrix onto the other subspace
        is_good_soln = false; lam1 = 0.1;
        while ~is_good_soln
            cvx_begin quiet
                variable betas(2*sum(N_o)+Ni,1) complex
                minimize(norm(D*I1*(xs+Q*betas),1)+lam1*norm(I2*(xs+Q*betas)))
                subject to
                norm(A_est*xs-b_nf+A_est*Q*betas) <= 1* fac1 * sig;
                norm(A_state*(xs+Q*betas))   <= fac_s;
            cvx_end
            if strcmp(cvx_status,'Solved')
                is_good_soln = true;
            elseif strcmp(cvx_status,'Inaccurate/Solved')
                fac1 = fac1*1.05;
                disp(['Inaccurate/solved, retrying with fac1 = ' num2str(fac1)]);
            else
                fac1 = fac1*1.05;
                disp(['Failed, retrying with fac1 = ' num2str(fac1)]);
            end
        end
        est_tf= xs + Q*betas;
    elseif mtd_id == 1
        disp('Solving with NM6');
        [U,S,V] = svd(A_est);
        
        %now implement morozov's principle
        fac = 1;
        p=1; xs = 0; convg = false;
        while ~convg
            xs = xs + V(:,p) * (U(:,p)'*b_nf)/S(p,p);
            res = norm(A_est*xs-b_nf);
            if res > fac*sig
                p = p+1;
            else
                convg = true;
            end
        end
        x1 = xs;
        
        Vds = V(:,1:p);
        
        [U1,S1,V1] = svd(A_state);
        
        %now implement morozov's principle
        fac_s = 0.1;
        p1=1; x2s = 0; convg = false;
        while ~convg
            x2s = x2s + V1(:,p1) * (U1(:,p1)'*b_state)/S1(p1,p1);
            res = norm(A_state*x2s-b_state);
            if res > fac_s
                p1 = p1+1;
            else
                convg = true;
            end
        end
        x_state = x2s;
        
        Vss = V1(:,1:p1);
        W = Vds'*Vss;
        [Uw,Sw,Vw] = svd(W);
        Wn = Vw(:,size(W,1)+1:end);
        P = orth(Vss*Wn);
        x2 = P*P'*x_state;
        
        fac1 = 0.7;
        Q = eye(size(A_est,2))-Vds*Vds'- P*P'; %the projector matrix onto the other subspace
        is_good_soln = false;
        while ~is_good_soln
            cvx_begin quiet
            variable x3(size(A_est,2),1) complex
            minimize(norm(D*(x1+x2+Q*x3),1))
            subject to
            norm(A_est*(x1+x2+Q*x3) - b_nf) <= 1* fac1 * sig;
            norm(A_state*(x1+x2+Q*x3) - b_state)   <= fac_s * norm(b_state);
            cvx_end
            if strcmp(cvx_status,'Solved')
                is_good_soln = true;
            elseif strcmp(cvx_status,'Inaccurate/Solved')
                fac1 = fac1*1.05;
                disp(['Inaccurate/solved, retrying with fac1 = ' num2str(fac1)]);
            else
                %              cvx_stat1{k} = {cvx_status};
                fac1 = fac1*1.05;
                disp(['Failed, retrying with fac1 = ' num2str(fac1)]);
            end
        end
        est_tf = x1+x2+Q*x3;
    elseif mtd_id == 2
        disp('Solving with ADMM');
        [U,S,V] = svd(A_est);
        %now implement morozov's principle
        fac = 1;
        p=1; xs = 0; convg = false;
        while ~convg
            xs = xs + V(:,p) * (U(:,p)'*b_nf)/S(p,p);
            res = norm(A_est*xs-b_nf);
            if res > fac*sig
                p = p+1;
            else
                convg = true;
            end
        end
        %now implement SOM and L1 minimization
        Vs = V(:,1:p);
        Q = eye(size(A_est,2))-Vs*Vs'; %the projector matrix onto the other subspace
        noise_var = 10^(-snr/10)* norm(b_nf);
        I1 = [eye(2*sum(N_o)),zeros(2*sum(N_o),Ni)]; I2 = [zeros(Ni,2*sum(N_o)),eye(Ni)];
        notconvg = true; i = 1; maxiter = 500;
        pnorm = 0.1; %value of p-norm
        lambda1 = 0.5; lambda2 = 1; nu = 0.1; %%%%%%%% Need to set %%%%%%%%%
        M = D(1:2*sum(N_o),1:2*sum(N_o));
        %keep reducing the lambda value till the mse goes below noise variance
        tic
        while notconvg
            lambda1 = lambda1 * 0.5;
            if i==1
                x_g = zeros(2*sum(N_o)+Ni,1); v = ones(2*sum(N_o)+Ni,1);
            else
                x_g = x_recov; %warm start
            end
            [x_recov,rmse] = incfield_admm(A_est,A_state,Q,x_g,xs,M,I1,I2,b_nf,pnorm,nu,lambda1,lambda2,maxiter);
            %figure; stem(x); hold;stem(x_recov); legend;grid;
            obj(i) = rmse; lambda_val(i) = lambda1;
            i = i+1;
            if rmse < 2*sqrt(3*noise_var)
                notconvg = false;
            end
        end
        toc
        est_tf = xs+Q*x_recov;
    end
    est_tf = est_tf./vnorm';
    terror = norm(est_tf(1:2*sum(N_o))-x_full_oni)/norm(x_full_oni);
    disp(['Tangential field error is ' num2str(terror)]);
    %generate prediction matrix
    
    if run_prd
        tic;
        ind = [];
        for i = 1:size(locs_g,1)
            r_g = locs_g(i,:);
            if in_air(r_g,invp,0.1*lambda)
                ind = [ind; i];
            end
        end
        locs_gp = locs_g(ind,:); %locations for prediction
        B_grid = gen_matrices(invp,locs_gp,false,inc_id); %Prediction matrix
        B_gridi = gen_incmtx(invp,locs_gp,Ni,false,inc_id);
        save(fname_e,'B_grid','B_gridi','locs_gp','dl_g','ind');
        disp(['Prediction matrix on grid to ' fname_g]);
        toc;
    else
        load(fname_e);
    end
    est_grid = [B_grid,B_gridi] * est_tf;
    gerror = norm(est_grid-sol_grid(ind)-sol_gridi(ind))...
        /norm(sol_grid(ind)+sol_gridi(ind));
    disp(['Reconstruction error is ' num2str(gerror)]);
    save(['result_' scheme '.mat'],'terror','gerror','est_grid','sol_grid','sol_gridi','ind','locs_gp','locs','invp');
    if plot_2d == true
        dl_g = lambda/20;
        extents = invp.corners(1).tr - invp.corners(1).bl;
        nx = floor(extents(1)/dl_g);
        ny = floor(extents(2)/dl_g);
        x_g = dl_g/2:dl_g:(nx-dl_g/2)*dl_g;
        y_g = fliplr(dl_g/2:dl_g:(ny-dl_g/2)*dl_g);
        [X_g,Y_g] = meshgrid(x_g,y_g);
        Epr = zeros(size(X_g));
        for i = 1:size(X_g,1)
            for j = 1:size(X_g,2)
                ind1 = find(locs_gp(:,1)==X_g(i,j) & locs_gp(:,2)==Y_g(i,j));
                if isempty(ind1)
                    Epr(i,j) = 0.01;
                else
                    Epr(i,j) = est_grid(ind1);
                end
            end
        end
        figure;
        imagesc(abs(Epr));
        toc;
    end
end

%% End of main code, Helper functions
function [b] = sampled_noise(y,snr)
%complex white Gaussian noise at the given per-sample SNR (dB); same as
%awgn(y(i),snr,'measured') without the Communications Toolbox
b=zeros(size(y));
for i=1:length(y)
    noise_power = abs(y(i))^2/10^(snr/10);
    b(i) = y(i) + sqrt(noise_power/2)*(randn + 1j*randn);
end
end
function legit = in_air(r,geom,excl)
%function returns true if point r is inside region 1
%optional exclusion zone around each object possible
legit = true;
%Check outer most boundary first
x_L = geom.corners(1).bl(1,1) + excl;
x_R = geom.corners(1).tr(1,1) - excl;
y_B = geom.corners(1).bl(1,2) + excl;
y_T = geom.corners(1).tr(1,2) - excl;
if r(1) < x_L || r(1) > x_R || r(2) < y_B || r(2) > y_T
    legit = false;
end
%Check object boundaries now
for j=2:length(geom.surf)
    if geom.shape(j,:)=='rect'
        tx_L = geom.corners(j).bl(1,1) - excl;
        tx_R = geom.corners(j).tr(1,1) + excl;
        ty_B = geom.corners(j).bl(1,2) - excl;
        ty_T = geom.corners(j).tr(1,2) + excl;
        if r(1)>=tx_L && r(1)<=tx_R && r(2)>=ty_B && r(2)<=ty_T
            legit = false;
        end
    else
        cen = geom.corners(j).bl;
        rad = geom.corners(j).tr(1,1);
        if norm(r-cen) < (rad + excl)
            legit = false;
        end
    end
end
% finally, check for source
tx_L = geom.source.bl(1,1) - excl;
tx_R = geom.source.tr(1,1) + excl;
ty_B = geom.source.bl(1,2) - excl;
ty_T = geom.source.tr(1,2) + excl;
if r(1)>=tx_L && r(1)<=tx_R && r(2)>=ty_B && r(2)<=ty_T
    legit = false;
end
end
function [A_i] = gen_incmtx(invp,locs,Ni,state,inc_id)
%Generate incident field submatrix for a given geometry and set of locations
%i.e. Grafs Addition theorem
nmeas = size(locs,1); %no of meas points
A_i = zeros(nmeas,1); %estimation matrix
eps = invp.eps(1); %called from inverse solver, only one eps of R1
kvec = 2*pi*sqrt(eps)/invp.lambda;
mu0 =  4*pi*10^(-7);
omega = kvec*3e8;
if inc_id == 0
    r0 = [invp.lambda/2+5*invp.lambda ...
        -3*invp.lambda/4+5*invp.lambda];
elseif inc_id == 1
    r0 = [5*invp.lambda 5*invp.lambda];
end

for i = 1:nmeas
    p = locs(i,:)-r0;
    for j = 1:Ni
            j1 = -(Ni-1)/2+j-1;
            A_i(i,j) = (-1)^(j1+1)*mu0*omega/4*besselh(j1,2,kvec*norm(p))*exp(1j*j1*angle(p(1)+1j*p(2)));
    end
end
if state
    A_i = -A_i; %in state eqn, inc field is subtracted
end
end
function [A_est,b_est] = gen_matrices(geom,locs,state,inc_id)
%Generate data/state matrix for a given geometry and set of locations
%i.e. Huygen's principle or extinction thm for Region 1 alone
%nreg = numel(geom.eps);%number of regions
%count number of variables, n(i) has the nvars of surf i-1
nvar = zeros(numel(geom.surf),1);
for i = 1:numel(geom.surf)
    nvar(i)= 2 * length(geom.surf(i).m(:,1));
end
N = sum(nvar);
nmeas = size(locs,1); %no of meas points
A_est = zeros(nmeas,N); %estimation matrix
[w,z] = glquadrule(2); %Use an n-pt quadrature rule
if numel(geom.eps) == 1
    eps = geom.eps(1); %called from inverse solver, only one eps of R1
else
    eps = geom.eps(2); %called from forward sover, one eps for each R
end
kvec = 2*pi*sqrt(eps)/geom.lambda;
for i=1:nmeas
    offset = 0;
    for j=1:length(geom.surf)
        for k = 1:length(geom.surf(j).m)
            [g,dg] = intg_green(locs(i,:),geom.surf(j).p(k,:),...
                geom.surf(j).q(k,:),geom.surf(j).n(k,:),kvec,1,w,z);
            A_est(i,offset+k) = -g;
            A_est(i,offset+k+length(geom.surf(j).m)) = +dg;
        end
        offset = offset + 2*k;
    end
end
b_est = [];
if state
    A_est = -A_est; %sign difference for state eqn terms
    b_est = zeros(nmeas,1);
    for i=1:nmeas
        b_est(i) = incfield(eps,geom.lambda,locs(i,:),inc_id);
    end
end

end
function [y] = incfield(epsr,lambda,r,inc_id)
%Get the inc field at a given point
k = 2*pi*sqrt(epsr)/lambda;
if inc_id == 0
    r0 = [lambda/2+5*lambda -3*lambda/4+5*lambda];
    y = besselh(0,2,k.*norm(r-r0)); %point cyl source
elseif inc_id == 1
    % J(x',y') = x'^2+y'^2-x'+y'   if norm((x',y'),inf) > lambda/8
    % J(x',y') = 0                else
    % Ng point gauss legendre in 2D
    r1 = r-[5*lambda,5*lambda];
    if (norm(r1,inf)<lambda/8)
        y = 0;
        return;
    end
    Ng = 6;
    omega = 2*pi*3e8/lambda;
    mu0 = 4*pi*1e-7;
    [x,w] = lgwt(Ng,-1,1); 
    [xx,yy,wtw] = lgwt2d(x,w); %Only within unit circle
    Ni = length(xx);
    y = 0; 
%     f1 = 2/(abs(besselh(0,2,20*pi))*625*pi^2);
    f1 = 1.5/(omega*mu0*abs(besselh(0,2,20*pi))*(lambda/8)^4);
    for j = 1:Ni
        p1 = (lambda/8)*[xx(j),yy(j)];
        y = y + f1*1j*omega*mu0*wtw(j)*((lambda/8)^2)...
            *(-1j/4*besselh(0,2,k*norm(r1-p1)))*(p1(1)^2+p1(2)^2-p1(1)+p1(2));
    end
end
end
function [xx,yy,ww] = lgwt2d(x,w)
    [X,Y]=meshgrid(x,x);
    [W1,W2] = meshgrid(w,w);
    ww = reshape(W1.*W2,[numel(W1),1]);
    xx = reshape(X,[numel(X),1]);     yy = reshape(Y,[numel(Y),1]);
end
function [g,dg] = intg_green(m,p,q,n,k,reg,w,z)
%This fn integrates green and gradgreen over a segment from p to q with
%normal vector n at a testing point m in a region 'reg' with wavevector k
%use the supplied quadrule weights (w) and knots (z)
%green = -j/4*besselh(0,2,k*rho),
%gradgreen = j*k/4*besselh(1,2,k*rho)\hat(rho).n, rho = r-m, r\in[p,q]
gam  = 0.5772;
h = norm(p-q);
if norm(m-(p+q)*0.5) < 1e-6 %the singularity terms
    g = -1i/4*h*(1-1i*2/pi*(log(k*h/4)+gam-1));
    if reg == 1
        dg = -0.5; %For Region 1 alone (inward normals)
    else
        dg =  0.5; %For all other regions.
    end
else
    %parameterize over seg by t \in [0,1] as q+(1-t)(p-q)
    green = @(t) -1i/4*besselh(0,2,k*norm(q+(1-t)*(p-q)-m));
    gradgreen = @(t) 1i*k/4*besselh(1,2,k*norm(q+(1-t)*(p-q)-m))* ...
        dot(n,q+(1-t)*(p-q)-m)/norm(q+(1-t)*(p-q)-m);
    %todo: move these fns out of here to speed up
    g = h*glquad(0,1,green,w,z);
    dg = h*glquad(0,1,gradgreen,w,z);
end
end
function [p,m,q,n] = gen_surfgrid(BL,TR,dl,vis,shp)
%generate surface grid for rect object: BL to TR
%generate surface grid for circ object: center at BL, radius TR(1)
%Generates midpt of coordinates and outward normals as size nx2
%try to get the best discretization close to given dl
%p is start of seg, m is middle, q is end, and n is normal
if shp == 'rect'
    ext = TR-BL;
    nx = floor(ext(1)/dl); ny = floor(ext(2)/dl);
    dl_x = ext(1)/nx; dl_y = ext(2)/ny;
    n = zeros(2*(nx+ny),2); %x,y coordinates of normal vector
    m = n; %x,y coordinates of midpoint
    p = n; q = n; %x,y coordinates of segment end points
    %set contour coordinates in clockwise manner
    %from BL to TL
    ind = 1:ny; offset = ny;
    m(ind,1) = 0; m(ind,2) = 0.5*dl_y:dl_y:(ny-0.5)*dl_y;
    p(ind,1) = 0; p(ind,2) = m(ind,2) - 0.5*dl_y;
    q(ind,1) = 0; q(ind,2) = m(ind,2) + 0.5*dl_y;
    n(ind,:) = repmat([-1,0],length(ind),1);
    %from TL to TR
    ind = offset+1:offset+nx; offset = offset+nx;
    m(ind,1) = 0.5*dl_x:dl_x:(nx-0.5)*dl_x; m(ind,2) = ny*dl_y;
    p(ind,1) = m(ind,1)-0.5*dl_x; p(ind,2) = m(ind,2);
    q(ind,1) = m(ind,1)+0.5*dl_x; q(ind,2) = m(ind,2);
    n(ind,:) = repmat([0,1],length(ind),1);
    %from TR to BR
    ind = offset+1:offset+ny; offset = offset+ny;
    m(ind,1) = nx*dl_x; m(ind,2) = (ny-0.5)*dl_y:-dl_y:0.5*dl_y;
    p(ind,1) = m(ind,1); p(ind,2) = m(ind,2)+0.5*dl_y;
    q(ind,1) = m(ind,1); q(ind,2) = m(ind,2)-0.5*dl_y;
    n(ind,:) = repmat([1,0],length(ind),1);
    %for BR to BL
    ind = offset+1:offset+nx;
    m(ind,1) = (nx-0.5)*dl_x:-dl_x:0.5*dl_x; m(ind,2) = 0;
    p(ind,1) = m(ind,1)+0.5*dl_x; p(ind,2) = m(ind,2);
    q(ind,1) = m(ind,1)-0.5*dl_x; q(ind,2) = m(ind,2);
    n(ind,:) = repmat([0,-1],length(ind),1);
    %move to desired origin
    m = m + BL; p = p + BL; q = q + BL;
elseif shp == 'circ'
    nc = floor(2*pi*TR(1)/dl);
    n = zeros(nc,2); %x,y coordinates of normal vector
    m = n; %x,y coordinates of midpoint
    p = n; q = n; %x,y coordinates of segment end points
    th = 0:2*pi/nc:2*pi-2*pi/nc;
    th_p = th-pi/nc; th_q = th+pi/nc;
    %         m(:,1) = BL(1) + TR(1)*cos(th'); m(:,2) = BL(2) + TR(1)*sin(th');
    p(:,1) = BL(1) + TR(1)*cos(th_p'); p(:,2) = BL(2) + TR(1)*sin(th_p');
    q(:,1) = BL(1) + TR(1)*cos(th_q'); q(:,2) = BL(2) + TR(1)*sin(th_q');
    m(:,1) = (p(:,1)+q(:,1))*0.5; m(:,2) = (p(:,2)+q(:,2))*0.5;
    tv = q-p; %tangential vectors
    for i = 1:size(tv,1)
        n(i,:) = [tv(i,2), -tv(i,1)]/norm(tv(i,:));
    end
end
%test visualization
if ~isempty(vis)
    quiver(vis,m(:,1),m(:,2),n(:,1),n(:,2)); %plot normals from mid-point
    scatter(vis,p(:,1),p(:,2),'rx');
    scatter(vis,q(:,1),q(:,2),'bo');
    axis equal;
end
end

