TO247_length = 16; % mm
TO247_height = 5; % mm
TO247_width = 21; % mm

TM_length = 16; % mm
TM_height = 2; % mm
TM_width = 21; % mm

chassis_radius = 21; % mm
chassis_depth = 50; % mm

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
    
    hi_probdef('millimeters','planar',1e-8,chassis_depth,30);

% Add material hi_addmaterial("materialname", kx, ky, qv, kt)
    %kx:Thermal conductivity in the x- or r-direction in unit of W/(m K)
    %ky:Thermal conductivity in the y- or z-direction in unit of W/(m K)
    %qv:Volume heat generation density in units of W/m3
    %kt:Volumetric heat capacity in units of MJ/(m3*K)

    hi_addmaterial('Silicon', 130, 130, Q_mosfet);
    hi_addmaterial('Aluminum', 239, 239, 0);
    hi_addmaterial('Ceramic', 25, 25, 0);
    
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

    hi_addboundprop('Top boundary', 2, 0, 0, 423, 0, 0);
    hi_addboundprop('Bot boundary', 2, 0, 0, 423, 2.3, 0);
    
    x1_temp = 0;
    y1_temp = 0;
    x2_temp = 2*chassis_radius;
    y2_temp = 0;
    x3_temp = 2*chassis_radius;
    y3_temp = chassis_radius;
    x4_temp = chassis_radius+TM_length/2;
    y4_temp = chassis_radius;
    x5_temp = chassis_radius+TM_length/2;
    y5_temp = chassis_radius+TM_height;
    x6_temp = chassis_radius+TM_length/2;
    y6_temp = chassis_radius+TM_height+TO247_height;
    x7_temp = chassis_radius-TM_length/2;
    y7_temp = chassis_radius+TM_height+TO247_height;
    x8_temp = chassis_radius-TM_length/2;
    y8_temp = chassis_radius+TM_height;
    x9_temp = chassis_radius-TM_length/2;
    y9_temp = chassis_radius;
    x10_temp = 0;
    y10_temp = chassis_radius;
    
    % 1
    hi_addnode(x1_temp,y1_temp);
    % 2
    hi_addnode(x2_temp,y2_temp);
    % 3
    hi_addnode(x3_temp,y3_temp);
    % 4
    hi_addnode(x4_temp,y4_temp);
    % 5
    hi_addnode(x5_temp,y5_temp);
    % 6
    hi_addnode(x6_temp,y6_temp);
    % 7
    hi_addnode(x7_temp,y7_temp);
    % 8
    hi_addnode(x8_temp,y8_temp);
    % 9
    hi_addnode(x9_temp,y9_temp);
    % 10
    hi_addnode(x10_temp,y10_temp);
    
% create geometry

    % chassis

    hi_addsegment(x1_temp,y1_temp,x2_temp,y2_temp);
    hi_addsegment(x2_temp,y2_temp,x3_temp,y3_temp);
    hi_addsegment(x3_temp,y3_temp,x4_temp,y4_temp);
    hi_addsegment(x4_temp,y4_temp,x9_temp,y9_temp);
    hi_addsegment(x9_temp,y9_temp,x10_temp,y10_temp);
    hi_addsegment(x10_temp,y10_temp,x1_temp,y1_temp);
    
    hi_addblocklabel(chassis_radius,chassis_radius/2);
    
    % apply the material 'Aluminum' to the chassis
    hi_selectlabel(chassis_radius,chassis_radius/2);
    % hi_setblockprop("blockname", automesh, meshsize, group)
    % automesh: 0 = mesher defers to mesh size constraint defined in meshsize, 1 = mesher
    % automatically chooses the mesh density.
    % meshsize: size constraint on the mesh in the block marked by this label.
    hi_setblockprop('Aluminum',1,0,0);
    hi_clearselected();
    
    % thermal interface, TM

    hi_addsegment(x4_temp,y4_temp,x5_temp,y5_temp);
    hi_addsegment(x5_temp,y5_temp,x8_temp,y8_temp);
    hi_addsegment(x8_temp,y8_temp,x9_temp,y9_temp);
    
    hi_addblocklabel(chassis_radius,chassis_radius + TM_height/2);
    
    % apply the material 'ceramic' to TM
    hi_selectlabel(chassis_radius,chassis_radius + TM_height/2);
    % hi_setblockprop("blockname", automesh, meshsize, group)
    % automesh: 0 = mesher defers to mesh size constraint defined in meshsize, 1 = mesher
    % automatically chooses the mesh density.
    % meshsize: size constraint on the mesh in the block marked by this label.
    hi_setblockprop('Ceramic',1,0,0);
    hi_clearselected();
    
    % MOSFET

    hi_addsegment(x5_temp,y5_temp,x6_temp,y6_temp);
    hi_addsegment(x6_temp,y6_temp,x7_temp,y7_temp);
    hi_addsegment(x7_temp,y7_temp,x8_temp,y8_temp);
    
    hi_addblocklabel(chassis_radius,chassis_radius + TM_height + TO247_height/2);
    
    % apply the material 'Silicon' to the MOSFET
    hi_selectlabel(chassis_radius,chassis_radius + TM_height + TO247_height/2);
    % hi_setblockprop("blockname", automesh, meshsize, group)
    % automesh: 0 = mesher defers to mesh size constraint defined in meshsize, 1 = mesher
    % automatically chooses the mesh density.
    % meshsize: size constraint on the mesh in the block marked by this label.
    hi_setblockprop('Silicon',1,0,0);
    hi_clearselected();
    
% Apply boundary condition

    hi_selectsegment(0,chassis_radius/2); 
    hi_selectsegment(2*chassis_radius,chassis_radius/2); 
    hi_selectsegment(chassis_radius,0);
    % hi setsegmentprop("propname", elementsize, automesh, hide, group, "inconductor")
        % Boundary property "propname"
        % Local element size along segment no greater than elementsize
        % automesh: 0 = mesher defers to the element constraint defined by elementsize, 1 = mesher
        % automatically chooses mesh size along the selected segments
        % hide: 0 = not hidden in post-processor, 1 == hidden in post processor
        % A member of group number group
        % A member of the conductor specified by the string "inconductor". If the segment is not
        % part of a conductor, this parameter can be specified as "<None>".
    hi_setsegmentprop('Bot boundary',0,1,0,0,'<None>');
    hi_clearselected;
    
    hi_selectsegment((chassis_radius - TM_length/2)/2,chassis_radius); 
    hi_selectsegment(2*chassis_radius - (chassis_radius - TM_length/2)/2,chassis_radius); 
    hi_selectsegment(chassis_radius,chassis_radius + TM_height + TO247_height);
    hi_selectsegment(chassis_radius - TM_length/2,chassis_radius + TM_height/2);
    hi_selectsegment(chassis_radius - TM_length/2,chassis_radius + TM_height + TO247_height/2);
    hi_selectsegment(chassis_radius + TM_length/2,chassis_radius + TM_height/2);
    hi_selectsegment(chassis_radius + TM_length/2,chassis_radius + TM_height + TO247_height/2);
    hi_setsegmentprop('Top boundary',0,1,0,0,'<None>');
    hi_clearselected;
    
% save file
    hi_saveas(['MOS_therm_' num2str(Pid.ID) '.feh']);
    
% analyze and load solution

    hi_analyze();
    hi_loadsolution();
    
% get max. case temperature of MOSFET
    
    Heat_sample_size = 1;
    MOSFET_x_sample = 0:1:TO247_length;
    MOSFET_T_array = zeros(1,size(MOSFET_x_sample,2));
    count_MOSFET_sample = 0;
    for i=1:1:size(MOSFET_x_sample,2)
        count_MOSFET_sample = count_MOSFET_sample + 1;
        T_temp = ho_getpointvalues(chassis_radius - ...
            TM_length/2+MOSFET_x_sample(1,i),chassis_radius + TM_height);
        MOSFET_T_array(1,count_MOSFET_sample) = T_temp(1);
    end
    
    % maximum case temperature of MOSFET
    T_MOSFET_max = max(MOSFET_T_array)- 273.15; % celsius degree
    
    closefemm;
    
    
    
    
    
    