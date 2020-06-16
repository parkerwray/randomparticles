
clc
addpath(genpath('/home/parkerwray/hypnos/Codes/Matlab Functions')); 
addpath(genpath('/home/parkerwray/hypnos/Codes/raacampbell-shadedErrorBar-9b19a7b')); 
%% PROCESS 3D RANDOM NANOPARTICLE FILMS
%{
This code processes simulations from Main_Random_3D_...
The first objective is to 

%}
lda = wavelengths;
[A, B, Ax, Bx] = convert_modes(sphere_coeffs_par, sphere_coeffs_perp);

%%

figure, 
plot(lda, mean(abs(A(:,:,1,1)),1))
hold on 
plot(lda, mean(abs(B(:,:,1,1)),1)) 
hold off



%%
figure, 
hold on 
plot(lda, mean(abs(A(:,:,1,1)),1))
plot(lda, mean(abs(B(:,:,1,1)),1))
plot(lda, mean(abs(A(:,:,1,2)),1))
plot(lda, mean(abs(B(:,:,1,2)),1))
plot(lda, mean(abs(A(:,:,2,1)),1))
plot(lda, mean(abs(B(:,:,2,1)),1))
plot(lda, mean(abs(A(:,:,2,2)),1))
plot(lda, mean(abs(B(:,:,2,2)),1))
plot(lda, mean(abs(A(:,:,2,3)),1))
plot(lda, mean(abs(B(:,:,2,3)),1))



%set(h_legend,'FontSize',34)
xlim([min(lda),max(lda)])
title(['Efficiency of a-Si nanoparticle, d = ', num2str(2*r1), 'nm'])
set(gca,'FontSize',34)








%convert_modes2(sphere_coeffs_par.')
% % % %% Process data of the entire particle cluster
% % % 
% % % GET_CLUSTER_DATA = 0;
% % % 
% % % if GET_CLUSTER_DATA  == 1
% % %     [cluster_file, FileName, PathName] = get_file();
% % %     load(cluster_file);
% % % end
% % % 
% % % lda = wavelengths;
% % % 
% % % [c_par, S_par] = ...
% % %     unpack_cluster(cluster_par);
% % % [c_perp, S_perp] = ...
% % %     unpack_cluster(cluster_perp);
% % % c = combine_cluster(c_perp, c_par);
% % % c = cluster_stats(c);
% % % 
% % % figure, 
% % % plot_avg_efficiencies(lda, ...
% % %     [c.qe_mean,c.qe_std],...
% % %     [c.qs_mean,c.qs_std],...
% % %     [c.qa_mean,c.qa_std]);
% % % title(strcat('3D Random film ', num2str(100*ff), 'FF Parallel'))
% % % 
% % % 
% % % %% Process efficiencies of each particle in cluster
% % % 
% % % % Get the absorption, scattering, and extinction cross-section of each
% % % % sphere in the MSTM simulation for each simulation. 
% % % s_par = unpack_sphere(sphere_result_par);
% % % s_perp = unpack_sphere(sphere_result_perp);
% % % s = combine_sphere(s_perp, s_par);
% % % s = sphere_stats(s);
% % % 
% % % i=5;
% % % figure, 
% % % plot_avg_efficiencies(lda, ...
% % %     [s{i}.qe_mean,s{i}.qe_std],...
% % %     [s{i}.qs_mean,s{i}.qs_std],...
% % %     [s{i}.qa_mean,s{i}.qa_std]);
% % % title(strcat('Avg. Particle Efficiency for ', num2str(100*ff),...
% % %     'FF for distribution ', num2str(i)))
% % % % Plot the efficiency result for each sphere for each simulation at a given
% % % % wavelength
% % % figure,
% % % hold on 
% % % for i = 1:size(s,2)
% % %     [value, idx] = min(abs(lda-390));
% % %     histogram(s{i}.qa(idx,:),round(size(s{1}.qa(idx,:),2)/4))
% % % end
% % % hold off


%% Process Modes of each particle in cluster










%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% SUB-FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_avg_efficiencies(lda, qe,qa,qs)
    hold on 
    errplot(lda, qe(:,1), qe(:,2))
    errplot(lda, qs(:,1), qs(:,2))
    errplot(lda, qa(:,1), qa(:,2))
    %errplot(lda, asym(:,1), asym(:,2))
    hold off
    ylabel('Efficiency')
    xlabel('Wavelength (nm)')
    legend('\sigma_{ext.}','\sigma_{sca.}','\sigma_{abs.}')
end


function errplot(x,y,error)
    s = shadedErrorBar(x, y, error,...
        'lineprops',{'LineWidth',6},'patchSaturation', 0.1);
    set(s.edge, 'LineWidth',0.5)
end


%% Functions for sphere efficiency processing

function [s] = unpack_sphere(dummy)
    
    for instance_idx = 1:size(dummy,2)
        for lda_idx = 1:size(dummy,1)
            dummy2 = dummy{lda_idx, instance_idx};
                
                s{instance_idx}.qe(lda_idx,:) = ...
                    cell2mat({dummy2.qe}); 
                s{instance_idx}.qa(lda_idx,:) = ...
                    cell2mat({dummy2.qa});                 
                s{instance_idx}.qs(lda_idx,:) = ...
                    cell2mat({dummy2.qs});                 
                s{instance_idx}.qav(lda_idx,:) = ...
                    cell2mat({dummy2.qav});                 
        end
    end
    
end

function s = combine_sphere(s1, s2)
    for idx = 1:size(s1,2)
        s{idx}.qe = [s1{idx}.qe,s2{idx}.qe];
        s{idx}.qa = [s1{idx}.qa,s2{idx}.qa];
        s{idx}.qs = [s1{idx}.qs,s2{idx}.qs];
        s{idx}.qav = [s1{idx}.qav,s2{idx}.qav];
    end
end

function s = sphere_stats(s)
    for idx = 1:size(s,2)
        s{idx}.qe_mean = mean(s{idx}.qe,2);
        s{idx}.qe_std = std(s{idx}.qe,[],2);
        s{idx}.qs_mean = mean(s{idx}.qs,2);
        s{idx}.qs_std = std(s{idx}.qs,[],2);    
        s{idx}.qa_mean = mean(s{idx}.qa,2);
        s{idx}.qa_std = std(s{idx}.qa,[],2);
        s{idx}.qav_mean = mean(s{idx}.qav,2);
        s{idx}.qav_std = std(s{idx}.qav,[],2);
    end
end


%% Functions for cluster level processing
function c = cluster_stats(c)
    c.qe_mean = mean(c.qe,2);
    c.qe_std = std(c.qe,[],2);
    c.qs_mean = mean(c.qs,2);
    c.qs_std = std(c.qs,[],2);    
    c.qa_mean = mean(c.qa,2);
    c.qa_std = std(c.qa,[],2);
    c.asym_mean = mean(c.asym,2);
    c.asym_std = std(c.asym,[],2);
end

function c = combine_cluster(c1, c2)

    c.qe = [c1.qe,c2.qe];
    c.qa = [c1.qa,c2.qa];
    c.qs = [c1.qs,c2.qs];
    c.asym = [c1.asym,c2.asym];

end

function [c,S] = unpack_cluster(cluster_results)



    for lda_idx = 1:size(cluster_results,1)
        for instance_idx = 1:size(cluster_results,2)
            % Grab all the coeffs for a specific lda, instance pair
            c.qe(lda_idx, instance_idx) = ...
                cluster_results{lda_idx, instance_idx}.qext; 
            c.qs(lda_idx, instance_idx) = ...
                cluster_results{lda_idx, instance_idx}.qsca;             
            c.qa(lda_idx, instance_idx) = ...
                cluster_results{lda_idx, instance_idx}.qabs;     
            c.asym(lda_idx, instance_idx) = ...
                cluster_results{lda_idx, instance_idx}.asym;  
            S{lda_idx, instance_idx} = ...
                cluster_results{lda_idx, instance_idx}.S; 
        end
    end
    
end

%% Functions for mode level processing


function [A, B, Ax, Bx] = unpack_modes(modes, n, polarization)
% Parallel polarization = 0
% Perpendicular polarization = 1
A = [];
B = [];
Ax = [];
Bx = [];

        orig = modes(1); % Grab the sphere at the origin data
        modes_a_par = orig.a_te{n};
        modes_b_par = orig.a_tm{n};

        %Sort Coeffs
        m = (-n:1:n);
        A_neg_m_n = [];
        B_neg_m_n = [];
        A_m_n = [];
        B_m_n = [];
        for idx_dummy = 1:length(m)
            if m(idx_dummy)<0
                A_neg_m_n = [A_neg_m_n,  modes_a_par(idx_dummy)];
                B_neg_m_n = [B_neg_m_n,  modes_b_par(idx_dummy)];
            elseif m(idx_dummy)>0
                A_m_n = [A_m_n, modes_a_par(idx_dummy)];
                B_m_n = [B_m_n, modes_b_par(idx_dummy)];
            elseif m(idx_dummy) == 0
                A_0_n = modes_a_par(idx_dummy);
                B_0_n = modes_b_par(idx_dummy);
            end
        end
        
        
        A = NaN(n+1,1);
        B = NaN(n+1,1);
        Ax = NaN(n+1,1);
        Bx = NaN(n+1,1);
        
        
        

        A = A_0_n;
        B = B_0_n;
        Ax = NaN;
        Bx = NaN;
        for m = 1:n
            dummy = ...
                A_m_n(m)+((-1).^(m+polarization)).*A_neg_m_n(m);

            A =...
                [A, dummy];

            dummy = ...
                B_m_n(m)+((-1).^(m+1+polarization)).*B_neg_m_n(m);

            B =...
                [B, 2.*dummy];

            dummy = ...
                A_m_n(m)+((-1).^(m+1+polarization)).*A_neg_m_n(m);

            Ax =...
                [Ax, dummy];

            dummy = ...
                B_m_n(m)+((-1).^(m+polarization)).*B_neg_m_n(m);

            Bx =...
                [Bx, 2.*dummy];                   


        end




end

function [A, B, Ax, Bx, Apar, Bpar, Axpar, Bxpar,...
            Aper, Bper, Axper, Bxper] = ...
                convert_modes(sphere_modes_par, sphere_modes_per)


if size(sphere_modes_par, 1) > 1
    [len_lda, len_sims] = size(sphere_modes_par);
else
    len_lda = length(sphere_modes_par);
    len_sims = 1;
end

Max_Order = 5;    

Apar= NaN(len_sims, len_lda, Max_Order, Max_Order+1); 
Bpar= NaN(len_sims, len_lda, Max_Order, Max_Order+1); 
Axpar= NaN(len_sims, len_lda, Max_Order, Max_Order+1); 
Bxpar= NaN(len_sims, len_lda, Max_Order, Max_Order+1); 
Aper= NaN(len_sims, len_lda, Max_Order, Max_Order+1); 
Bper= NaN(len_sims, len_lda, Max_Order, Max_Order+1); 
Axper= NaN(len_sims, len_lda, Max_Order, Max_Order+1); 
Bxper= NaN(len_sims, len_lda, Max_Order, Max_Order+1); 

% Get Parallel Modes
    for idx_sims = 1:len_sims
        for idx_lda = 1:len_lda
            Nmax = sphere_modes_par{idx_lda,idx_sims}.Nsphere_order;
            for idx_n = 1:Nmax
            polarization = 0; % Flag for parallel is 0
            [A,B,Ax,Bx] = ...
                unpack_modes(sphere_modes_par{idx_lda,idx_sims},idx_n,...
                                polarization);    
                for idx_m = 1:length(A)     
                    if isempty(A(idx_m))
                        Apar(idx_sims, idx_lda, idx_n, idx_m) = NaN;
                    else
                        Apar(idx_sims, idx_lda, idx_n, idx_m) = A(idx_m);
                    end  
                    if isempty(B(idx_m))
                        Bpar(idx_sims, idx_lda, idx_n, idx_m) = NaN;
                    else
                        Bpar(idx_sims, idx_lda, idx_n, idx_m) = B(idx_m);
                    end                      
                    if isempty(Ax(idx_m))
                        Axpar(idx_sims, idx_lda, idx_n, idx_m) = NaN;
                    else
                        Axpar(idx_sims, idx_lda, idx_n, idx_m) = Ax(idx_m);
                    end  
                    if isempty(Bx(idx_m))
                        Bxpar(idx_sims, idx_lda, idx_n, idx_m) = NaN;
                    else
                        Bxpar(idx_sims, idx_lda, idx_n, idx_m) = Bx(idx_m);
                    end                        
                    
                end
            
            end
        end
    end
% Get Perpendicular Modes
    for idx_sims = 1:len_sims
        for idx_lda = 1:len_lda
            Nmax = sphere_modes_par{idx_lda,idx_sims}.Nsphere_order;
            for idx_n = 1:Nmax
            polarization = 1; % Flag for perpendicular is 1
            [A,B,Ax,Bx] = ...
                unpack_modes(sphere_modes_per{idx_lda,idx_sims},idx_n,...
                                polarization);
                for idx_m = 1:length(A)     
                    if isempty(A(idx_m))
                        Aper(idx_sims, idx_lda, idx_n, idx_m) = NaN;
                    else
                        Aper(idx_sims, idx_lda, idx_n, idx_m) = A(idx_m);
                    end  
                    if isempty(B(idx_m))
                        Bper(idx_sims, idx_lda, idx_n, idx_m) = NaN;
                    else
                        Bper(idx_sims, idx_lda, idx_n, idx_m) = B(idx_m);
                    end                      
                    if isempty(Ax(idx_m))
                        Axper(idx_sims, idx_lda, idx_n, idx_m) = NaN;
                    else
                        Axper(idx_sims, idx_lda, idx_n, idx_m) = Ax(idx_m);
                    end  
                    if isempty(Bx(idx_m))
                        Bxper(idx_sims, idx_lda, idx_n, idx_m) = NaN;
                    else
                        Bxper(idx_sims, idx_lda, idx_n, idx_m) = Bx(idx_m);
                    end                        
                    
                end
            end
        end
    end
    
    A = cat(1, Apar, Aper);
    B = cat(1, Bpar, Bper);
    Ax = cat(1, Axpar, Axper);
    Bx = cat(1, Bxpar, Bxper);
    

%     A1par = cell2mat({Apar{1, :, 1}}.');
%     B1par = cell2mat({Bpar{1, :, 1}}.');
%     A2par = cell2mat({Apar{1, :, 2}}.');
%     B2par = cell2mat({Bpar{1, :, 2}}.');
%     A3par = cell2mat({Apar{1, :, 3}}.');
%     B3par = cell2mat({Bpar{1, :, 3}}.');
%     
%     figure, 
%     hold on 
%     plot(abs(A1par).^2)
%     plot(abs(B1par).^2)
%     plot(abs(A2par).^2)
%     plot(abs(B2par).^2)
%     plot(abs(A3par).^2)
%     plot(abs(B3par).^2)    
%     hold off
%     legend('a01 par', 'a11 par', 'b01 par','b11 par',...
%         'a02 par', 'a12 par','a22 par','b02 par','b12 par', 'b22 par',...
%         'a03 par', 'a13 par','a23 par','a33 par','b03 par','b13 par','b23 par', 'b33 par');
%     Ax1par = cell2mat({Axpar{1, :, 1}}.');
%     Bx1par = cell2mat({Bxpar{1, :, 1}}.');
%     Ax2par = cell2mat({Axpar{1, :, 2}}.');
%     Bx2par = cell2mat({Bxpar{1, :, 2}}.');
%     
%     figure, 
%     hold on 
%     plot(abs(Ax1par).^2)
%     plot(abs(Bx1par).^2)
%     plot(abs(Ax2par).^2)
%     plot(abs(Bx2par).^2)
%     hold off
%     legend('ax11 par', 'bx11 par',...
%         'ax12 par','ax22 par','bx12 par','bx22 par')  
% 
%     A1per = cell2mat({Aper{1, :, 1}}.');
%     B1per = cell2mat({Bper{1, :, 1}}.');
%     A2per = cell2mat({Aper{1, :, 2}}.');
%     B2per = cell2mat({Bper{1, :, 2}}.');
% 
%     figure, 
%     hold on 
%     plot(abs(A1per).^2)
%     plot(abs(B1per).^2)
%     plot(abs(A2per).^2)
%     plot(abs(B2per).^2)
%     hold off
%     legend('a01 per', 'a11 per', 'b01 per','b11 per',...
%         'a02 per', 'a12 per','a22 per','b02 per','b12 per','b22 per')
%     
%     Ax1per = cell2mat({Axper{1, :, 1}}.');
%     Bx1per = cell2mat({Bxper{1, :, 1}}.');
%     Ax2per = cell2mat({Axper{1, :, 2}}.');
%     Bx2per = cell2mat({Bxper{1, :, 2}}.');
%     
%     figure, 
%     hold on 
%     plot(abs(Ax1per).^2)
%     plot(abs(Bx1per).^2)
%     plot(abs(Ax2per).^2)
%     plot(abs(Bx2per).^2)
%     hold off
%     legend( 'ax11 per', 'bx11 per',...
%          'ax12 per','ax22 per','bx12 per','bx22 per')      

end




function [anm1, anm2, bnm1, bnm2, x,...
    qext_par, qsca_par, qabs_par,...
    qext_per, qsca_per, qabs_per]= convert_modes2(sphere_coeffs)

%[-22,-12,02,12,22]
%[22,12,02,-12,-22]
%n = 2
% take [n+1,end] of both arrays
% add acording to m

dummy_te = {};
dummy_tm = {};

    for lda_idx = 1:size(sphere_coeffs,1)
        for instance_idx = 1:size(sphere_coeffs,2)
            % Grab all the coeffs for a specific lda, instance pair
            dummy_te = sphere_coeffs{lda_idx, instance_idx}.a_te;
            dummy_tm = sphere_coeffs{lda_idx, instance_idx}.a_tm;   
            x(lda_idx, instance_idx) = sphere_coeffs{lda_idx, instance_idx}.ka;
            qsca_par(lda_idx, instance_idx) = sphere_coeffs{lda_idx, instance_idx}.Qsca_par; 
            qabs_par(lda_idx, instance_idx) = sphere_coeffs{lda_idx, instance_idx}.Qabs_par;             
            qext_par(lda_idx, instance_idx) = sphere_coeffs{lda_idx, instance_idx}.Qext_par;             
            qsca_per(lda_idx, instance_idx) = sphere_coeffs{lda_idx, instance_idx}.Qsca_per; 
            qabs_per(lda_idx, instance_idx) = sphere_coeffs{lda_idx, instance_idx}.Qabs_per;             
            qext_per(lda_idx, instance_idx) = sphere_coeffs{lda_idx, instance_idx}.Qext_per; 
            
            
            
            for n = 1:length(dummy_te)
                if n == 1
                   % [-11,01,11]
                   dummy_tm_n = dummy_te{n};
                   dummy_te_n = dummy_tm{n};

                   a111 = (dummy_tm_n(3)-dummy_tm_n(1));
                   a011 = dummy_tm_n(2);

                   a112 = 2.*(dummy_te_n(3)-dummy_te_n(1));
                   a012 = dummy_te_n(2); 

                   b111 = (dummy_tm_n(3)+dummy_tm_n(1));
                   b011 = dummy_tm_n(2);

                   b112 = 2.*(dummy_te_n(3)+dummy_te_n(1));
                   b012 = dummy_te_n(2);               
                end

                if n == 2
                   %[-22,-12,02,12,22]                
                   dummy_tm_n = dummy_te{n};
                   dummy_te_n = dummy_tm{n};

                   a221 = (dummy_tm_n(5)+dummy_tm_n(1));
                   a121 = (dummy_tm_n(4)-dummy_tm_n(2));              
                   a021 = dummy_tm_n(3);

                   a222 = 2.*(dummy_te_n(5)+dummy_te_n(1));
                   a122 = 2.*(dummy_te_n(4)-dummy_te_n(2));              
                   a022 = dummy_te_n(3);               


                   b221 = (dummy_tm_n(5)-dummy_tm_n(1));
                   b121 = (dummy_tm_n(4)+dummy_tm_n(2));              
                   b021 = dummy_tm_n(3);

                   b222 = 2.*(dummy_te_n(5)-dummy_te_n(1));
                   b122 = 2.*(dummy_te_n(4)+dummy_te_n(2));              
                   b022 = dummy_te_n(3);              
                end            

                if n == 3
                   %[-33,-23,-13,03,13,23,33]                
                   dummy_tm_n = dummy_te{n};
                   dummy_te_n = dummy_tm{n};

                   a331 = (dummy_tm_n(7)-dummy_tm_n(1));
                   a231 = (dummy_tm_n(6)+dummy_tm_n(2));
                   a131 = (dummy_tm_n(5)-dummy_tm_n(1));              
                   a031 = dummy_tm_n(4);

                   a332 = 2.*(dummy_te_n(7)-dummy_te_n(1));
                   a232 = 2.*(dummy_te_n(6)+dummy_te_n(2));
                   a132 = 2.*(dummy_te_n(5)-dummy_te_n(1));              
                   a032 = dummy_te_n(4);              

                   b331 = (dummy_tm_n(7)+dummy_tm_n(1));
                   b231 = (dummy_tm_n(6)-dummy_tm_n(2));
                   b131 = (dummy_tm_n(5)+dummy_tm_n(1));              
                   b031 = dummy_tm_n(4);

                   b332 = 2.*(dummy_te_n(7)+dummy_te_n(1));
                   b232 = 2.*(dummy_te_n(6)-dummy_te_n(2));
                   b132 = 2.*(dummy_te_n(5)+dummy_te_n(1));              
                   b032 = dummy_te_n(4);                             
                end               
            end

            anm1(:,lda_idx,instance_idx) = [a011,a111,a021,a121,a221,a031,a131,a231,a331]; 
            anm2(:,lda_idx,instance_idx) = [a012,a112,a022,a122,a222,a032,a132,a232,a332]; 
            bnm1(:,lda_idx,instance_idx) = [b011,b111,b021,b121,b221,b031,b131,b231,b331]; 
            bnm2(:,lda_idx,instance_idx) = [b012,b112,b022,b122,b222,b032,b132,b232,b332];         

        end
    end


%     figure, 
%     hold on 
%     plot(abs(anm2(1,:,1)).^2)
%     plot(abs(anm2(2,:,1)).^2)
%     plot(abs(bnm2(1,:,1)).^2)
%     plot(abs(bnm2(2,:,1)).^2)
%     hold off
    figure, 
    hold on 
    plot(abs(anm1(1,:,1)).^2)
    plot(abs(anm1(2,:,1)).^2)
    plot(abs(bnm1(1,:,1)).^2)
    plot(abs(bnm1(2,:,1)).^2)
    plot(abs(anm2(1,:,1)).^2)
    plot(abs(anm2(2,:,1)).^2)
    plot(abs(bnm2(1,:,1)).^2)
    plot(abs(bnm2(2,:,1)).^2)
    
    hold off    
    legend('a011', 'a111', 'b011','b111', 'a012', 'a112', 'b012','b112')
    title('correct')
end




