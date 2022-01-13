% open femm with main window hiden
    openfemm(1);    
    %hand = HandleToFEMM;
% create a new preprocessor document: heat flow problem

    newdocument(2);
    
% define the problem
    % hi_probdef(units,type,precision,(depth),(minangle))

    % units: e 'inches', 'millimeters', 'centimeters', 'mils', 'meters', 
    % and 'micrometers'

    % type:  'planar' - 2-D planar problem, 'axi' - axisymmetric problem

    % precision: maixmum of RMS of the residual

    % depth: the depth of the problem in the into-the-page direction for
    % 2-D planar problems

    % minanle:  the minimum angle constraint sent to the mesh generator
    
    hi_probdef('millimeters','planar',1e-8,E_width_FEMM*N_stack_FEMM,30);

% Add material hi_addmaterial("materialname", kx, ky, qv, kt)
    %kx:Thermal conductivity in the x- or r-direction in unit of W/(m K)
    %ky:Thermal conductivity in the y- or z-direction in unit of W/(m K)
    %qv:Volume heat generation density in units of W/m3
    %kt:Volumetric heat capacity in units of MJ/(m3*K)

    hi_addmaterial('Iron', 4, 4, Q_core);
    hi_addmaterial('Copper', 386, 386, Q_winding);
    hi_addmaterial('FR4', 0.25, 0.25, 0);
    
% Add boundary conditions hi_addboundprop("boundpropname", BdryFormat, Tset, qs, Tinf, h, beta) 
    % To obtain a “Heat Flux” type boundary condition, set qs to be the heat flux density and
    % BdryFormat to 1. Set all other parameters to zero.
    
    % To obtain a convection boundary condition, set h to the desired heat transfer coefficient
    % and Tinf to the desired external temperature and set BdryFormat to 2.
    
    % For a Radiation boundary condition, set beta equal to the desired emissivity and Tinf
    % to the desired external temperature and set BdryFormat to 3.
    
    % For a “Periodic” boundary condition, set BdryFormat to 4 and set all other parameters
    % to zero.
    
    % For an “Anti-Perodic” boundary condition, set BdryFormat to 5 set all other parameters
    % to zero.

    hi_addboundprop('Outer boundary', 2, 0, 0, 423, 2.3, 0);
    
% create geometry

    % core
    hi_addnode(0,0);
    hi_addnode(E_length_O_FEMM / 2,0);
    hi_addnode(E_length_O_FEMM / 2,2 * E_height_O_FEMM);
    hi_addnode(0,2 * E_height_O_FEMM);
    
    hi_addnode(E_mid_FEMM / 2,E_height_O_FEMM - E_height_I_FEMM);
    hi_addnode(E_mid_FEMM / 2 + Window_width,E_height_O_FEMM - E_height_I_FEMM);
    hi_addnode(E_mid_FEMM / 2 + Window_width,2 * E_height_O_FEMM - E_height_I_FEMM);
    hi_addnode(E_mid_FEMM / 2,2 * E_height_O_FEMM - E_height_I_FEMM);
    
    hi_addsegment(0,0,E_length_O_FEMM / 2,0);
    hi_addsegment(E_length_O_FEMM / 2,0,E_length_O_FEMM / 2,2 * E_height_O_FEMM);
    hi_addsegment(E_length_O_FEMM / 2,2 * E_height_O_FEMM,0,2 * E_height_O_FEMM);
    hi_addsegment(0,2 * E_height_O_FEMM,0,0);
    
    hi_addsegment(E_mid_FEMM / 2,E_height_O_FEMM - E_height_I_FEMM,...
        E_mid_FEMM / 2 + Window_width,E_height_O_FEMM - E_height_I_FEMM);
    hi_addsegment(E_mid_FEMM / 2 + Window_width,E_height_O_FEMM - E_height_I_FEMM,...
        E_mid_FEMM / 2 + Window_width,2 * E_height_O_FEMM - E_height_I_FEMM);
    hi_addsegment(E_mid_FEMM / 2 + Window_width,2 * E_height_O_FEMM - E_height_I_FEMM,...
        E_mid_FEMM / 2,2 * E_height_O_FEMM - E_height_I_FEMM);
    hi_addsegment(E_mid_FEMM / 2,2 * E_height_O_FEMM - E_height_I_FEMM,...
        E_mid_FEMM / 2,E_height_O_FEMM - E_height_I_FEMM);
    
    hi_addblocklabel(E_mid_FEMM / 4,(E_height_O_FEMM - E_height_I_FEMM)/2);
    
    % apply the material 'Iron' to the core
    hi_selectlabel(E_mid_FEMM / 4,(E_height_O_FEMM - E_height_I_FEMM)/2);
    % hi_setblockprop("blockname", automesh, meshsize, group)
    % automesh: 0 = mesher defers to mesh size constraint defined in meshsize, 1 = mesher
    % automatically chooses the mesh density.
    % meshsize: size constraint on the mesh in the block marked by this label.
    hi_setblockprop('Iron',1,0,0);
    hi_clearselected();
    
    % winding
    for j = 1:1:N_B
        
        % primary winding
        if mod(j,2)
            for i = 1:1:N_P_L_FEMM 
                
                x1_temp = E_mid_FEMM / 2 + D_Horiz_P + (i - 1) * (winding_P_width + Winding_P_dist);
                y1_temp = E_height_O_FEMM - E_height_I_FEMM + D_vertical - winding_th_FEMM + (j - 1) * (D_layer_layer_FEMM + D_vertical);
                x2_temp = E_mid_FEMM / 2 + D_Horiz_P + winding_P_width + (i - 1) * (winding_P_width + Winding_P_dist);
                y2_temp = E_height_O_FEMM - E_height_I_FEMM + D_vertical + (j - 1) * (D_layer_layer_FEMM + D_vertical);
                
                hi_addnode(x1_temp,y1_temp);
                hi_addnode(x2_temp,y1_temp);
                hi_addnode(x2_temp,y2_temp);
                hi_addnode(x1_temp,y2_temp);
                
                hi_addsegment(x1_temp,y1_temp,x2_temp,y1_temp);
                hi_addsegment(x2_temp,y1_temp,x2_temp,y2_temp);
                hi_addsegment(x2_temp,y2_temp,x1_temp,y2_temp);
                hi_addsegment(x1_temp,y2_temp,x1_temp,y1_temp);

                hi_addblocklabel(E_mid_FEMM / 2 + D_Horiz_P + winding_P_width / 2 + (i - 1) * (winding_P_width + Winding_P_dist),...
                    E_height_O_FEMM - E_height_I_FEMM + D_vertical - winding_th_FEMM / 2 + (j - 1) * (D_layer_layer_FEMM + D_vertical));

                hi_selectlabel(E_mid_FEMM / 2 + D_Horiz_P + winding_P_width / 2 + (i - 1) * (winding_P_width + Winding_P_dist),...
                    E_height_O_FEMM - E_height_I_FEMM + D_vertical - winding_th_FEMM / 2 + (j - 1) * (D_layer_layer_FEMM + D_vertical));
                hi_setblockprop('Copper',1,0,0);
                hi_clearselected();

            end
        else
            for i = 1:1:N_P_L_FEMM 

                x1_temp = E_mid_FEMM / 2 + D_Horiz_P + (i - 1) * (winding_P_width + Winding_P_dist);
                y1_temp = E_height_O_FEMM - E_height_I_FEMM + D_vertical - winding_th_FEMM + (j - 1) * (D_layer_layer_FEMM + D_vertical) + D_layer_layer_FEMM + winding_th_FEMM;
                x2_temp = E_mid_FEMM / 2 + D_Horiz_P + winding_P_width + (i - 1) * (winding_P_width + Winding_P_dist);
                y2_temp = E_height_O_FEMM - E_height_I_FEMM + D_vertical + (j - 1) * (D_layer_layer_FEMM + D_vertical) + D_layer_layer_FEMM + winding_th_FEMM;
                
                hi_addnode(x1_temp,y1_temp);
                hi_addnode(x2_temp,y1_temp);
                hi_addnode(x2_temp,y2_temp);
                hi_addnode(x1_temp,y2_temp);
                
                hi_addsegment(x1_temp,y1_temp,x2_temp,y1_temp);
                hi_addsegment(x2_temp,y1_temp,x2_temp,y2_temp);
                hi_addsegment(x2_temp,y2_temp,x1_temp,y2_temp);
                hi_addsegment(x1_temp,y2_temp,x1_temp,y1_temp);

                hi_addblocklabel(E_mid_FEMM / 2 + D_Horiz_P + winding_P_width / 2 + (i - 1) * (winding_P_width + Winding_P_dist),...
                    E_height_O_FEMM - E_height_I_FEMM + D_vertical - winding_th_FEMM / 2 + (j - 1) * (D_layer_layer_FEMM + D_vertical) + D_layer_layer_FEMM + winding_th_FEMM);

                hi_selectlabel(E_mid_FEMM / 2 + D_Horiz_P + winding_P_width / 2 + (i - 1) * (winding_P_width + Winding_P_dist),...
                    E_height_O_FEMM - E_height_I_FEMM + D_vertical - winding_th_FEMM / 2 + (j - 1) * (D_layer_layer_FEMM + D_vertical) + D_layer_layer_FEMM + winding_th_FEMM);
                hi_setblockprop('Copper',1,0,0);
                hi_clearselected();

            end
        end

        % secondary winding
        if mod(j,2)
            for i = 1:1:N_S_L_FEMM 

                x1_temp = E_mid_FEMM / 2 + D_Horiz_S + (i - 1) * (winding_S_width + Winding_S_dist);
                y1_temp = E_height_O_FEMM - E_height_I_FEMM + D_vertical + D_layer_layer_FEMM + (j - 1) * (D_layer_layer_FEMM + D_vertical);
                x2_temp = E_mid_FEMM / 2 + D_Horiz_S + winding_S_width + (i - 1) * (winding_S_width + Winding_S_dist);
                y2_temp = E_height_O_FEMM - E_height_I_FEMM + D_vertical + D_layer_layer_FEMM + winding_th_FEMM + (j - 1) * (D_layer_layer_FEMM + D_vertical);
                
                hi_addnode(x1_temp,y1_temp);
                hi_addnode(x2_temp,y1_temp);
                hi_addnode(x2_temp,y2_temp);
                hi_addnode(x1_temp,y2_temp);
                
                hi_addsegment(x1_temp,y1_temp,x2_temp,y1_temp);
                hi_addsegment(x2_temp,y1_temp,x2_temp,y2_temp);
                hi_addsegment(x2_temp,y2_temp,x1_temp,y2_temp);
                hi_addsegment(x1_temp,y2_temp,x1_temp,y1_temp);

                hi_addblocklabel(E_mid_FEMM / 2 + D_Horiz_S + winding_S_width / 2 + (i - 1) * (winding_S_width + Winding_S_dist),...
                    E_height_O_FEMM - E_height_I_FEMM + D_vertical + D_layer_layer_FEMM + winding_th_FEMM / 2 + (j - 1) * (D_layer_layer_FEMM + D_vertical));

                hi_selectlabel(E_mid_FEMM / 2 + D_Horiz_S + winding_S_width / 2 + (i - 1) * (winding_S_width + Winding_S_dist),...
                    E_height_O_FEMM - E_height_I_FEMM + D_vertical + D_layer_layer_FEMM + winding_th_FEMM / 2 + (j - 1) * (D_layer_layer_FEMM + D_vertical));
                hi_setblockprop('Copper',1,0,0);
                hi_clearselected();

            end
        else
            for i = 1:1:N_S_L_FEMM 

                x1_temp = E_mid_FEMM / 2 + D_Horiz_S + (i - 1) * (winding_S_width + Winding_S_dist);
                y1_temp = E_height_O_FEMM - E_height_I_FEMM + j * (D_layer_layer_FEMM + D_vertical) - D_layer_layer_FEMM - winding_th_FEMM;
                x2_temp = E_mid_FEMM / 2 + D_Horiz_S + winding_S_width + (i - 1) * (winding_S_width + Winding_S_dist);
                y2_temp = E_height_O_FEMM - E_height_I_FEMM + winding_th_FEMM + j * (D_layer_layer_FEMM + D_vertical) - D_layer_layer_FEMM - winding_th_FEMM;
                
                hi_addnode(x1_temp,y1_temp);
                hi_addnode(x2_temp,y1_temp);
                hi_addnode(x2_temp,y2_temp);
                hi_addnode(x1_temp,y2_temp);
                
                hi_addsegment(x1_temp,y1_temp,x2_temp,y1_temp);
                hi_addsegment(x2_temp,y1_temp,x2_temp,y2_temp);
                hi_addsegment(x2_temp,y2_temp,x1_temp,y2_temp);
                hi_addsegment(x1_temp,y2_temp,x1_temp,y1_temp);

                hi_addblocklabel(E_mid_FEMM / 2 + D_Horiz_S + winding_S_width / 2 + (i - 1) * (winding_S_width + Winding_S_dist),...
                    E_height_O_FEMM - E_height_I_FEMM + winding_th_FEMM / 2 + j * (D_layer_layer_FEMM + D_vertical) - D_layer_layer_FEMM - winding_th_FEMM);

                hi_selectlabel(E_mid_FEMM / 2 + D_Horiz_S + winding_S_width / 2 + (i - 1) * (winding_S_width + Winding_S_dist),...
                    E_height_O_FEMM - E_height_I_FEMM + winding_th_FEMM / 2 + j * (D_layer_layer_FEMM + D_vertical) - D_layer_layer_FEMM - winding_th_FEMM);
                hi_setblockprop('Copper',1,0,0);
                hi_clearselected();

            end
        end
        
    
    end
    % FR4
                
    hi_addblocklabel(E_mid_FEMM / 2 + D_Horiz_P/2,...
        E_height_O_FEMM - E_height_I_FEMM + (D_vertical - winding_th_FEMM)/2);

    hi_selectlabel(E_mid_FEMM / 2 + D_Horiz_P/2,...
        E_height_O_FEMM - E_height_I_FEMM + (D_vertical - winding_th_FEMM)/2);
    hi_setblockprop('FR4',1,0,0);
    hi_clearselected();
    
% Apply boundary condition

    hi_selectsegment(0,E_height_O_FEMM); 
    hi_selectsegment(E_length_O_FEMM / 2,E_height_O_FEMM); 
    hi_selectsegment(E_length_O_FEMM / 4,0);
    hi_selectsegment(E_length_O_FEMM / 4,2 * E_height_O_FEMM);
    % hi setsegmentprop("propname", elementsize, automesh, hide, group, "inconductor")
        % Boundary property "propname"
        % Local element size along segment no greater than elementsize
        % automesh: 0 = mesher defers to the element constraint defined by elementsize, 1 = mesher
        % automatically chooses mesh size along the selected segments
        % hide: 0 = not hidden in post-processor, 1 == hidden in post processor
        % A member of group number group
        % A member of the conductor specified by the string "inconductor". If the segment is not
        % part of a conductor, this parameter can be specified as "<None>".
    hi_setsegmentprop('Outer boundary',0,1,0,0,'<None>');
    hi_clearselected;
    
% save file
    hi_saveas(['Tr_therm_' num2str(Pid.ID) '.feh']);
    
% analyze and load solution

    hi_analyze();
    hi_loadsolution();
    
% get max. case temperature of core
    
    Heat_sample_size = 1;
    Core_x_sample = 0:1:E_length_O_FEMM / 2;
    Core_y_sample = 0:1:2 * E_height_O_FEMM;
    Core_T_array = zeros(1,size(Core_x_sample,2)+size(Core_y_sample,2));
    count_core_sample = 0;
    for i=1:1:size(Core_x_sample,2)
        count_core_sample = count_core_sample + 1;
        Core_T_array(1,count_core_sample) = ho_getpointvalues(Core_x_sample(1,i),2 * E_height_O_FEMM);
    end
    for i=1:1:size(Core_y_sample,2)
        count_core_sample = count_core_sample + 1;
        Core_T_array(1,count_core_sample) = ho_getpointvalues(E_length_O_FEMM / 2,Core_y_sample(1,i));
    end
    
    % maximum case temperature of core
    T_core_max = max(Core_T_array);
    
    
    
    
    
    