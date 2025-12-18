clc
clear all
close all

% % global Constants

icount = 0;
aa =[0,3,4,8,11,12,16,19,20,24];
for ii = 1:10
    ig2 = aa(ii);
    a = sqrt(aa(ii));
    ia=int32(a);
    ic = 0;
    for igx = -ia:ia
        igx2 = igx*igx;
        for jgy=-ia:ia
            jgy2=jgy*jgy;
            igxy2=igx2+jgy2;
            for kgz=-ia:ia
                kgz2=kgz*kgz;
                igxyz2=igxy2+kgz2;
                if (igxyz2 == ig2)
                    ic=ic+1;
                    igx_vec(ic+icount)=igx;
                    igy_vec(ic+icount)=jgy;
                    igz_vec(ic+icount)=kgz;
                end
            end
        end
    end
    icount = ic + icount;
end

function E = calculate_Eigen(lattice_constant, v_sym, v_asym, path,igx_vec,igy_vec,igz_vec)
Q = 1.6e-19;  % charge in Colombs
h = 1.054571817e-34; % reduced planks constant(j)
mo = 9.11e-31; % kg is the free electron rest mass
M_e_eff = 1 * mo; %kg


%calculatiing off diagonal elements
S = (lattice_constant/8).*[1 1 1];
for row = 1:length(igx_vec)
    for col = 1:length(igy_vec)
        if row == col
            continue
        end
        diff = [double(igx_vec(row)-igx_vec(col)) double(igy_vec(row)-igy_vec(col)) double(igz_vec(row)-igz_vec(col))];
        mag_diff_sq = diff(1)^2 + diff(2)^2 + diff(3)^2;
        scalar_multi = diff(1)*S(1) + diff(2)*S(2) + diff(3)*S(3);
        temp2 = 1i*v_asym(mag_diff_sq)*sin(2*pi*scalar_multi/lattice_constant);
        temp1 = v_sym(mag_diff_sq)*cos(2*pi*scalar_multi/lattice_constant);
        H(row,col) =  temp1 + temp2;
    end
end

H_complete = zeros(length(H),length(H),length(path));

%calculating diagonal element
for width = 1:length(path)
    H_complete (:,:,width) = H;
    for row = 1:length(igx_vec)
        %addvec = [path(width,1)+double((igx_vec(row))*(2*pi/lattice_constant)) path(width,2)+double((igy_vec(row))*(2*pi/lattice_constant)) path(width,3)+double((igz_vec(row))*(2*pi/lattice_constant))];
        addvec = path(width,:) + double([igx_vec(row) igy_vec(row) igz_vec(row)]) *(2*pi/lattice_constant);
        mag_add_sq = addvec(1)^2  + addvec(2)^2 + addvec(3)^2;
        H_complete (row,row,width) = (h^2/(2*M_e_eff)) * mag_add_sq*(1/Q); %did not divide by Q as reduced planks constant is already in eV
    end
    eigenvalues(width,:)=(sort(eig(H_complete(:,:,width))))';
end

E = eigenvalues(:,1:10);
end

v_s = zeros(9,137);
v_a = zeros(9,137);
ry = 13.6; %ry to eV

%for Si @ 1
v_s(1,3) = -0.210 *ry;
v_s(1,8) = 0.040 *ry;
v_s(1,11) = 0.080 *ry;

%for Ge @ 2

v_s(2,3) = -0.230 *ry;
v_s(2,8) = 0.010 *ry;
v_s(2,11) = 0.060 *ry;

%for GaAs @ 3

v_s(3,3) = -0.230 *ry;
v_s(3,8) = 0.010 *ry;
v_s(3,11) = 0.060 *ry;

v_a(3,3) = 0.070*ry;
v_a(3,4) = 0.050*ry;
v_a(3,11) = 0.01*ry;

p = 0:0.2:1;

%difference between the v_s value GaAs and AlAs
diff_vs3 = (-0.221 + (0.230))/5;
diff_vs8 = (0.025 - (0.010))/5;
diff_vs11 = (0.070 - (0.060))/5;

diff_va3 = (0.080 - (0.070))/5;
diff_va4 = (0.050 - (0.050))/5;
diff_va11 = (-0.004 - (0.010))/5;

l_c = [5.43e-10 5.66e-10 5.64e-10];

for increment = 1:6
    v_s(increment+3, 3) = (-0.230 + (diff_vs3*(increment-1)))*ry;
    v_s(increment+3, 8) = (0.010 + (diff_vs8*(increment-1)))*ry;
    v_s(increment+3, 11) = (0.060 + (diff_vs11*(increment-1)))*ry;

    v_a(increment+3, 3) = (0.070 + (diff_va3*(increment-1)))*ry;
    v_a(increment+3, 4) = (0.050 + (diff_va4*(increment-1)))*ry;
    v_a(increment+3, 11) = (0.010 + (diff_va11*(increment-1)))*ry;

    l_c(increment+3) = (5.6533 + (0.0078*p(increment)))*(10^(-10));
    %defining the points values
end


for increment = 1:9
L = [pi/l_c(increment) pi/l_c(increment) pi/l_c(increment)];
gamma = [0 0 0];
X = [0 2*pi/l_c(increment) 0];
K = [(3/2)*(pi/l_c(increment)) (3/2)*(pi/l_c(increment)) 0];
ns = 20;

%finding ranges
drange1 = (gamma-L)/ns;
drange2 = (X-gamma)/ns;
drange3 = (K-X)/ns;
drange4 = (gamma-K)/ns;

%initial path allocations
path = zeros(ns*4,3);
path(1,:) = L;
path(ns+1,:) = gamma;
path(2*ns+1,:) = X;
path(3*ns+1,:) = K;

for count = 2:ns
    path(count,:) = path(count-1,:) + drange1;
    path(ns+count,:) = path(ns+count-1,:) + drange2;
    path((2*ns)+count,:) = path(2*ns+count-1,:) + drange3;
    path((3*ns)+count,:) = path(3*ns+count-1,:) + drange4;
end

    reduced_eigenvalues(increment,:,:) = calculate_Eigen(l_c(increment),v_s(increment,:),v_a(increment,:),path,igx_vec,igy_vec,igz_vec);

end


%minvalcb=zeros(1,9);
for increment = 1:9
    temperoryvar1 = squeeze(reduced_eigenvalues(increment,:,5));
   [minvalcb(increment), minindexcb(increment)] = min(temperoryvar1);
    temperoryvar2 = squeeze(reduced_eigenvalues(increment,:,4));
    [maxvalvb(increment),maxindexvb(increment)] = max(temperoryvar2);
    Band_gap(increment) = minvalcb(increment)-maxvalvb(increment); 
end

 figure (1)
 subplot(2,2,1)
 
 plot(squeeze(reduced_eigenvalues(1,:,:)),'HandleVisibility','off')
 hold on
 plot(minindexcb(1),minvalcb(1),'squarer','LineWidth',2,'DisplayName', "C.B Minima")
 plot(maxindexvb(1),maxvalvb(1),'^b','LineWidth',2,'DisplayName', "V.B Maxima")
 text(minindexcb(1)-0.2,minvalcb(1)+1,"\downarrow "+minvalcb(1)+"eV",'FontSize',9,'FontWeight','bold')
 text(maxindexvb(1)-0.2,maxvalvb(1)-1,"\uparrow "+maxvalvb(1)+"eV",'FontSize',9,'FontWeight','bold')
 plot(NaN,NaN, 'LineStyle', 'none','DisplayName', "E_{band gap} = " + Band_gap(1) + "eV" )
 legend show
grid on
hold off
 title ('Energy band structure of Si')
  xticks([1,11,21,31,40])
 xticklabels({'L','\gamma','X','K','\gamma'})
  xlabel('K')
 ylabel ('Energy (eV)')
 subplot(2,2,2)

 plot(squeeze(reduced_eigenvalues(2,:,:)),'HandleVisibility','off')
  hold on
 plot(minindexcb(2),minvalcb(2),'squarer','LineWidth',2,'DisplayName', "C.B Minima")
 plot(maxindexvb(2),maxvalvb(2),'^b','LineWidth',2,'DisplayName', "V.B Maxima")
 text(minindexcb(2)-0.2,minvalcb(2)+1,"\downarrow "+minvalcb(2)+"eV",'FontSize',9,'FontWeight','bold')
 text(maxindexvb(2)-0.2,maxvalvb(2)-1,"\uparrow "+maxvalvb(2)+"eV",'FontSize',9,'FontWeight','bold')
 plot(NaN,NaN, 'LineStyle', 'none','DisplayName', "E_{band gap}  = " + Band_gap(2) + "eV" )
 legend show
grid on
hold off
 title ('Energy band structure of Ge')
  xticks([1,11,21,31,40])
 xticklabels({'L','\gamma','X','K','\gamma'})
 xlabel('K')
 ylabel ('Energy (eV)')
 
 %figure(2)
 subplot(2,2,3)

 plot(squeeze(reduced_eigenvalues(3,:,:)),'HandleVisibility','off')
 hold on
 plot(minindexcb(3),minvalcb(3),'squarer','LineWidth',2,'DisplayName', "C.B Minima")
 plot(maxindexvb(3),maxvalvb(3),'^b','LineWidth',2,'DisplayName', "V.B Maxima")
 text(minindexcb(3)-0.2,minvalcb(3)+1,"\downarrow "+minvalcb(3)+"eV",'FontSize',9,'FontWeight','bold')
 text(maxindexvb(3)-0.2,maxvalvb(3)-1,"\uparrow "+maxvalvb(3)+"eV",'FontSize',9,'FontWeight','bold')
 plot(NaN,NaN, 'LineStyle', 'none','DisplayName', "E_{B.G} = " + Band_gap(3) + "eV" )
 legend show
grid on
hold off
title ('Energy band structure of GaAs')
 xticks([1,11,21,31,40])
 xticklabels({'L','\gamma','X','K','\gamma'})
 xlabel('K')
 ylabel ('Energy (eV)')
 
 subplot(2,2,4)
 
 plot(squeeze(reduced_eigenvalues(4,:,:)),'HandleVisibility','off')
  hold on
 plot(minindexcb(4),minvalcb(4),'squarer','LineWidth',2,'DisplayName', "C.B Minima")
 plot(maxindexvb(4),maxvalvb(4),'^b','LineWidth',2,'DisplayName', "V.B Maxima")
 text(minindexcb(4)-0.2,minvalcb(4)+1,"\downarrow "+minvalcb(4)+"eV",'FontSize',9,'FontWeight','bold')
 text(maxindexvb(4)-0.2,maxvalvb(4)-1,"\uparrow "+maxvalvb(4)+"eV",'FontSize',9,'FontWeight','bold')
 plot(NaN,NaN, 'LineStyle', 'none','DisplayName', "E_{B.G}  = " + Band_gap(4) + "eV" )
 legend show
grid on
hold off
 title ('Energy band structure of Al_{0}Ga_{1}As')
  xticks([1,11,21,31,40])
 xticklabels({'L','\gamma','X','K','\gamma'})
  xlabel('K')
 ylabel ('Energy (eV)')

  figure(3)
 subplot(2,2,1)
 plot(squeeze(reduced_eigenvalues(5,:,:)),'HandleVisibility','off')
  hold on
 plot(minindexcb(5),minvalcb(5),'squarer','LineWidth',2,'DisplayName', "C.B Minima")
 plot(maxindexvb(5),maxvalvb(5),'^b','LineWidth',2,'DisplayName', "V.B Maxima")
 text(minindexcb(5)-0.2,minvalcb(5)+1,"\downarrow "+minvalcb(5)+"eV",'FontSize',9,'FontWeight','bold')
 text(maxindexvb(5)-0.2,maxvalvb(5)-1,"\uparrow "+maxvalvb(5)+"eV",'FontSize',9,'FontWeight','bold')
 plot(NaN,NaN, 'LineStyle', 'none','DisplayName', "E_{B.G} = " + Band_gap(5) + "eV" )
 legend show
grid on
hold off
 title ('Energy band structure of Al_{0.2}Ga_{0.8}As')
  xticks([1,11,21,31,40])
 xticklabels({'L','\gamma','X','K','\gamma'})
 xlabel('K')
 ylabel ('Energy (eV)')
 subplot(2,2,2)
 plot(squeeze(reduced_eigenvalues(6,:,:)),'HandleVisibility','off')
  hold on
 plot(minindexcb(6),minvalcb(6),'squarer','LineWidth',2,'DisplayName', "C.B Minima")
 plot(maxindexvb(6),maxvalvb(6),'^b','LineWidth',2,'DisplayName', "V.B Maxima")
 text(minindexcb(6)-0.2,minvalcb(6)+1,"\downarrow "+minvalcb(1)+"eV",'FontSize',9,'FontWeight','bold')
 text(maxindexvb(6)-0.2,maxvalvb(6)-1,"\uparrow "+maxvalvb(1)+"eV",'FontSize',9,'FontWeight','bold')
 plot(NaN,NaN, 'LineStyle', 'none','DisplayName', "E_{B.G} = " + Band_gap(6) + "eV" )
 legend show
grid on
hold off
 title ('Energy band structure of Al_{0.4}Ga_{0.6}As')
  xticks([1,11,21,31,40])
 xticklabels({'L','\gamma','X','K','\gamma'})
 xlabel('K')
 ylabel ('Energy (eV)')

   %figure(4)
 subplot(2,2,3)
 plot(squeeze(reduced_eigenvalues(7,:,:)),'HandleVisibility','off')
 hold on
 plot(minindexcb(7),minvalcb(7),'squarer','LineWidth',2,'DisplayName', "C.B Minima")
 plot(maxindexvb(7),maxvalvb(7),'^b','LineWidth',2,'DisplayName', "V.B Maxima")
 text(minindexcb(7)-0.2,minvalcb(7)+1,"\downarrow "+minvalcb(7)+"eV",'FontSize',9,'FontWeight','bold')
 text(maxindexvb(7)-0.2,maxvalvb(7)-1,"\uparrow "+maxvalvb(7)+"eV",'FontSize',9,'FontWeight','bold')
 plot(NaN,NaN, 'LineStyle', 'none','DisplayName', "E_{B.G} = " + Band_gap(7) + "eV" )
 legend show
grid on
hold off 
 title ('Energy band structure of Al_{0.6}Ga_{0.4}As')
  xticks([1,11,21,31,40])
 xticklabels({'L','\gamma','X','K','\gamma'})
  xlabel('K')
 ylabel ('Energy (eV)')
 subplot(2,2,4)
 plot(squeeze(reduced_eigenvalues(8,:,:)),'HandleVisibility','off')
 hold on
 plot(minindexcb(8),minvalcb(8),'squarer','LineWidth',2,'DisplayName', "C.B Minima")
 plot(maxindexvb(8),maxvalvb(8),'^b','LineWidth',2,'DisplayName', "V.B Maxima")
 text(minindexcb(8)-0.2,minvalcb(8)+1,"\downarrow "+minvalcb(8)+"eV",'FontSize',9,'FontWeight','bold')
 text(maxindexvb(8)-0.2,maxvalvb(8)-1,"\uparrow "+maxvalvb(8)+"eV",'FontSize',9,'FontWeight','bold')
 plot(NaN,NaN, 'LineStyle', 'none','DisplayName', "E_{B.G} = " + Band_gap(8) + "eV" )
 legend show
grid on
hold off 
 title ('Energy band structure of Al_{0.8}Ga_{0.2}As')
  xticks([1,11,21,31,40])
 xticklabels({'L','\gamma','X','K','\gamma'})
 xlabel('K')
 ylabel ('Energy (eV)')

   figure(5)
 plot(squeeze(reduced_eigenvalues(9,:,:)),'HandleVisibility','off')
 hold on
 plot(minindexcb(9),minvalcb(9),'squarer','LineWidth',2,'DisplayName', "C.B Minima")
 plot(maxindexvb(9),maxvalvb(9),'^b','LineWidth',2,'DisplayName', "V.B Maxima")
 text(minindexcb(9)-0.2,minvalcb(9)+1,"\downarrow "+minvalcb(9)+"eV",'FontSize',9,'FontWeight','bold')
 text(maxindexvb(9)-0.2,maxvalvb(9)-1,"\uparrow "+maxvalvb(9)+"eV",'FontSize',9,'FontWeight','bold')
 plot(NaN,NaN, 'LineStyle', 'none','DisplayName', "E_{B.G} = " + Band_gap(9) + "eV" )
 legend show
grid on
hold off 
 title ('Energy band structure of Al_{1}Ga_{0}As')
  xticks([1,11,21,31,40])
 xticklabels({'L','\gamma','X','K','\gamma'})
 xlabel('K')
 ylabel ('Energy (eV)')

 figure(6)
 plot(p,Band_gap(4:length(Band_gap)))
 xlabel('Value of x')
 ylabel('Energy Band gap (eV)')
 title('Bandgap (E_{g}) variation of Al_{x}Ga_{1-x}As with alloy concentration (x)')