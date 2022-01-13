% open femm with main window hiden
    openfemm(1);    
    %hand = HandleToFEMM;
% create a new preprocessor document: electrostatics problem

    newdocument(1);

% define the problem
    % ei_probdef(units,type,precision,depth,minangle)

    % units: e 'inches', 'millimeters', 'centimeters', 'mils', 'meters', 
    % and 'micrometers'

    % type:  'planar' - 2-D planar problem, 'axi' - axisymmetric problem

    % precision: maixmum of RMS of the residual

    % depth: the depth of the problem in the into-the-page direction for
    % 2-D planar problems

    % minable:  the minimum angle constraint sent to the mesh generator

    ei_probdef('millimeters','planar',1e-8,E_width_FEMM*N_stack_FEMM,30);
    
% Add material
    % ei_addmaterial(’matname’, ex, ey, qv) adds a new material with called ’matname’
    % with the material properties:
    % ex Relative permittivity in the x- or r-direction.
    % ey Relative permittivity in the y- or z-direction.
    % qv Volume charge density in units of C/m3
    
    ei_addmaterial('FR4', 4, 4, 0);
    
% Define boundary condition
    % ei_addboundprop(’boundname’, Vs, qs, c0, c1, BdryFormat) adds a new boundary property with name ’boundname’
    % For a “Fixed Voltage” type boundary condition, set the Vs parameter to the desired voltage
    % and all other parameters to zero.
    % To obtain a “Mixed” type boundary condition, set C1 and C0 as required and BdryFormat to
    % 1. Set all other parameters to zero.
    % To obtain a prescribes surface charge density, set qs to the desired charge density in C/m2
    % and set BdryFormat to 2.
    % For a “Periodic” boundary condition, set BdryFormat to 3 and set all other parameters to zero.
    % For an “Anti-Perodic” boundary condition, set BdryFormat to 4 set all other parameters to zero.
    ei_addboundprop('ground', 0, 0, 0, 0, 0);

% create geometry

    % core
    
    % inner profile
    ei_drawrectangle([E_mid_FEMM / 2 E_height_O_FEMM - E_height_I_FEMM; E_length_I_FEMM / 2 E_height_O_FEMM + E_height_I_FEMM]);
    
    ei_addconductorprop('core', 0,0,1);
    
    ei_addblocklabel(E_length_O_FEMM / 4,E_height_O_FEMM - E_height_I_FEMM);
    ei_addblocklabel(E_length_O_FEMM / 4,E_height_O_FEMM + E_height_I_FEMM);
    ei_addblocklabel(E_mid_FEMM / 2,E_height_O_FEMM);
    ei_addblocklabel(E_length_I_FEMM / 2,E_height_O_FEMM);
    
    ei_selectlabel(E_length_O_FEMM / 4,E_height_O_FEMM - E_height_I_FEMM);
    ei_selectlabel(E_length_O_FEMM / 4,E_height_O_FEMM + E_height_I_FEMM);
    ei_selectlabel(E_mid_FEMM / 2,E_height_O_FEMM);
    ei_selectlabel(E_length_I_FEMM / 2,E_height_O_FEMM);
    
    % ei_setsegmentprop(’name’, elmsize, automesh, hide, group, ’inconductor’)
    % Set the select segments to have:
    % Boundary property ’name’
    % Local element size along segment no greater than elmsize
    % automesh: 0 = mesher defers to the element constraint defined by elementsize, 1 = mesher
    % automatically chooses mesh size along the selected segments
    % hide: 0 = not hidden in post-processor, 1 == hidden in post processor
    % A member of group number group
    % A member of the conductor specified by the string ’inconductor’. If the segment is not
    % part of a conductor, this parameter can be specified as ’<None>’.
    ei_setsegmentprop('ground', 0, 1, 1, 0, 'core');
    ei_clearselected;
    
    % winding
    for j = 1:1:N_B
        
        % primary winding
        if mod(j,2)
            for i = 1:1:N_P_L_FEMM 

                ei_drawrectangle([E_mid_FEMM / 2 + D_Horiz_P + (i - 1) * (winding_P_width + Winding_P_dist) ...
                    E_height_O_FEMM - E_height_I_FEMM + D_vertical - winding_th_FEMM + (j - 1) * (D_layer_layer_FEMM + D_vertical); ...
                    E_mid_FEMM / 2 + D_Horiz_P + winding_P_width + (i - 1) * (winding_P_width + Winding_P_dist) ...
                    E_height_O_FEMM - E_height_I_FEMM + D_vertical + (j - 1) * (D_layer_layer_FEMM + D_vertical)]);

                ei_addblocklabel(E_mid_FEMM / 2 + D_Horiz_P + winding_P_width / 2 + (i - 1) * (winding_P_width + Winding_P_dist),...
                    E_height_O_FEMM - E_height_I_FEMM + D_vertical - winding_th_FEMM / 2 + (j - 1) * (D_layer_layer_FEMM + D_vertical));
                
                % ei_addconductorprop(’conductorname’, Vc, qc, conductortype)
                % adds a new conductor property with name ’conductorname’ with either a prescribed 
                % voltage or a prescribed total charge. Set the unused property to zero. The conductortype parameter is 0
                % for prescribed charge and 1 for prescribed voltage.
                ei_addconductorprop(['coil_B_' num2str(j) '_P_' num2str(i)], I_P_FEMM,0,1);

                ei_selectlabel(E_mid_FEMM / 2 + D_Horiz_P + winding_P_width / 2 + (i - 1) * (winding_P_width + Winding_P_dist),...
                    E_height_O_FEMM - E_height_I_FEMM + D_vertical - winding_th_FEMM / 2 + (j - 1) * (D_layer_layer_FEMM + D_vertical));
                ei_setblockprop('Copper', 0, 1, ['coil_B_' num2str(j) '_P_' num2str(i)], 0, 0, 0);
                ei_clearselected;

            end
        else
            for i = 1:1:N_P_L_FEMM 

                mi_drawrectangle([E_mid_FEMM / 2 + D_Horiz_P + (i - 1) * (winding_P_width + Winding_P_dist) ...
                    E_height_O_FEMM - E_height_I_FEMM + D_vertical - winding_th_FEMM + (j - 1) * (D_layer_layer_FEMM + D_vertical) + D_layer_layer_FEMM + winding_th_FEMM; ...
                    E_mid_FEMM / 2 + D_Horiz_P + winding_P_width + (i - 1) * (winding_P_width + Winding_P_dist) ...
                    E_height_O_FEMM - E_height_I_FEMM + D_vertical + (j - 1) * (D_layer_layer_FEMM + D_vertical) + D_layer_layer_FEMM + winding_th_FEMM]);

                mi_addblocklabel(E_mid_FEMM / 2 + D_Horiz_P + winding_P_width / 2 + (i - 1) * (winding_P_width + Winding_P_dist),...
                    E_height_O_FEMM - E_height_I_FEMM + D_vertical - winding_th_FEMM / 2 + (j - 1) * (D_layer_layer_FEMM + D_vertical) + D_layer_layer_FEMM + winding_th_FEMM);

                mi_addcircprop(['coil_B_' num2str(j) '_P_' num2str(i)], I_P_FEMM, 1);

                mi_selectlabel(E_mid_FEMM / 2 + D_Horiz_P + winding_P_width / 2 + (i - 1) * (winding_P_width + Winding_P_dist),...
                    E_height_O_FEMM - E_height_I_FEMM + D_vertical - winding_th_FEMM / 2 + (j - 1) * (D_layer_layer_FEMM + D_vertical) + D_layer_layer_FEMM + winding_th_FEMM);
                mi_setblockprop('Copper', 0, 1, ['coil_B_' num2str(j) '_P_' num2str(i)], 0, 0, 0);
                mi_clearselected;

            end
        end

        % secondary winding
        if mod(j,2)
            for i = 1:1:N_S_L_FEMM 

                mi_drawrectangle([E_mid_FEMM / 2 + D_Horiz_S + (i - 1) * (winding_S_width + Winding_S_dist) ...
                    E_height_O_FEMM - E_height_I_FEMM + D_vertical + D_layer_layer_FEMM + (j - 1) * (D_layer_layer_FEMM + D_vertical); ...
                    E_mid_FEMM / 2 + D_Horiz_S + winding_S_width + (i - 1) * (winding_S_width + Winding_S_dist) ...
                    E_height_O_FEMM - E_height_I_FEMM + D_vertical + D_layer_layer_FEMM + winding_th_FEMM + (j - 1) * (D_layer_layer_FEMM + D_vertical)]);

                mi_addblocklabel(E_mid_FEMM / 2 + D_Horiz_S + winding_S_width / 2 + (i - 1) * (winding_S_width + Winding_S_dist),...
                    E_height_O_FEMM - E_height_I_FEMM + D_vertical + D_layer_layer_FEMM + winding_th_FEMM / 2 + (j - 1) * (D_layer_layer_FEMM + D_vertical));

                mi_addcircprop(['coil_B_' num2str(j) '_S_' num2str(i)], - I_S, 1);

                mi_selectlabel(E_mid_FEMM / 2 + D_Horiz_S + winding_S_width / 2 + (i - 1) * (winding_S_width + Winding_S_dist),...
                    E_height_O_FEMM - E_height_I_FEMM + D_vertical + D_layer_layer_FEMM + winding_th_FEMM / 2 + (j - 1) * (D_layer_layer_FEMM + D_vertical));
                mi_setblockprop('Copper', 0, 1, ['coil_B_' num2str(j) '_S_' num2str(i)], 0, 0, 0);
                mi_clearselected;

            end
        else
            for i = 1:1:N_S_L_FEMM 

                mi_drawrectangle([E_mid_FEMM / 2 + D_Horiz_S + (i - 1) * (winding_S_width + Winding_S_dist) ...
                    E_height_O_FEMM - E_height_I_FEMM + j * (D_layer_layer_FEMM + D_vertical) - D_layer_layer_FEMM - winding_th_FEMM; ...
                    E_mid_FEMM / 2 + D_Horiz_S + winding_S_width + (i - 1) * (winding_S_width + Winding_S_dist) ...
                    E_height_O_FEMM - E_height_I_FEMM + winding_th_FEMM + j * (D_layer_layer_FEMM + D_vertical) - D_layer_layer_FEMM - winding_th_FEMM]);

                mi_addblocklabel(E_mid_FEMM / 2 + D_Horiz_S + winding_S_width / 2 + (i - 1) * (winding_S_width + Winding_S_dist),...
                    E_height_O_FEMM - E_height_I_FEMM + winding_th_FEMM / 2 + j * (D_layer_layer_FEMM + D_vertical) - D_layer_layer_FEMM - winding_th_FEMM);

                mi_addcircprop(['coil_B_' num2str(j) '_S_' num2str(i)], - I_S, 1);

                mi_selectlabel(E_mid_FEMM / 2 + D_Horiz_S + winding_S_width / 2 + (i - 1) * (winding_S_width + Winding_S_dist),...
                    E_height_O_FEMM - E_height_I_FEMM + winding_th_FEMM / 2 + j * (D_layer_layer_FEMM + D_vertical) - D_layer_layer_FEMM - winding_th_FEMM);
                mi_setblockprop('Copper', 0, 1, ['coil_B_' num2str(j) '_S_' num2str(i)], 0, 0, 0);
                mi_clearselected;

            end
        end
        
    
    end
    
    % air space, as boundary
    mi_drawrectangle([-margin -margin; E_length_O_FEMM / 2 + margin ...
        2 * E_height_O_FEMM + margin]);
    
% Add block labels

    % core
    
    
    % air window
    mi_addblocklabel(E_length_O_FEMM / 4,E_height_O_FEMM - E_height_I_FEMM + (D_vertical - winding_th_FEMM) / 2);
    
    % air space, as boundary
    mi_addblocklabel(- margin / 2,- margin / 2);
    
    % block labels for coil have been added when creating geometry
    
% Add circuit property, have done when creating geometry
    
% Apply the materials to the appropriate block label

    mi_selectlabel(- margin / 2,- margin / 2);
    mi_setblockprop('Air', 0, 1, '<None>', 0, 0, 0);
    mi_clearselected;
    
    mi_selectlabel(E_length_O_FEMM / 4,E_height_O_FEMM - E_height_I_FEMM + (D_vertical - winding_th_FEMM) / 2);
    mi_setblockprop('Air', 0, 1, '<None>', 0, 0, 0);
    mi_clearselected;
    
    
    
% save file
    mi_saveas(['Tr_FCT_' num2str(Pid.ID) '.fem']);
    
% analyze and load solution

    mi_analyze;
    mi_loadsolution;
    
%     
% get current, voltage, and flux linkage of coils

    primary_coils = zeros(N_P_L_FEMM * N_B,4);
    secondary_coils = zeros(N_S_L_FEMM * N_B,4);

    for j = 1:1:N_B
        
        % primary winding
        for i = 1:1:N_P_L_FEMM 

            vals = mo_getcircuitproperties(['coil_B_' num2str(j) '_P_' num2str(i)]);
            primary_coils(i + N_P_L_FEMM * (j - 1),1) = vals(1);
            primary_coils(i + N_P_L_FEMM * (j - 1),2) = vals(2);
            primary_coils(i + N_P_L_FEMM * (j - 1),3) = vals(3);
            primary_coils(i + N_P_L_FEMM * (j - 1),4) = real(vals(2) / vals(1));

        end

        % secondary winding
        for i = 1:1:N_S_L_FEMM 

            vals = mo_getcircuitproperties(['coil_B_' num2str(j) '_S_' num2str(i)]);
            secondary_coils(i + N_S_L_FEMM * (j - 1),1) = vals(1);
            secondary_coils(i + N_S_L_FEMM * (j - 1),2) = vals(2);
            secondary_coils(i + N_S_L_FEMM * (j - 1),3) = vals(3);
            secondary_coils(i + N_S_L_FEMM * (j - 1),4) = real(vals(2) / vals(1));

        end
    
    end
        
    R_pri = sum(primary_coils(:,4));
    R_sec = sum(secondary_coils(:,4));
    Loss = I_P_FEMM * I_P_FEMM * R_pri + I_S * I_S * R_sec;
%     
    closefemm;