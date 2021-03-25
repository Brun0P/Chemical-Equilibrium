classdef Compound
    %Compound
    % 
    
    properties
        NameCompound  % Name of the molecule -> "water" or "butane"
        Atoms         % Array 4x1 with nº of carbon, hydrogen, oxygen and nitrogen
        HoT           % Standard Enthalpy J/mol
        GoT           % Standard Gibbs Free Energy J/mol
        Cp            % Array with Cp coefficients
        Tc            % Critic temperature
        Pc            % Critic pressure
        Vc            % Critic volume
        Zc            % Z...
        Omega         % Accentric factor
    end
    
    methods
        
        function obj = Compound(NameCompound)
            %Compound Constructor
            %
            
            [obj.Atoms, obj.HoT, obj.GoT, obj.Cp, obj.Tc, obj.Pc, obj.Vc, obj.Zc, obj.Omega] = obj.AssignProperties(NameCompound);
            obj.NameCompound = NameCompound;
        end
        
        function [Atoms, HoT, GoT, Cp, Tc, Pc, Vc, Zc, Omega] = AssignProperties(obj, Name)
            %AssignAtoms
            %   Given the compound name, it is assigned an array 4x1 with
            %   nº of carbon, hydrogen, oxygen and nitrogen, respectively.
            %   Also, assigns the properties of each molecule.
            
            % Cp equation
            % Cp/R = A + BT + CT^2 + DT^-2 + ET^3
            
            switch Name % We can change this to SMILES string format
                case 'Carbon monoxide'
                    Atoms = [1,0,1,0];
                    Prop = [-110.53	 -137.16   3.376 0.557 0 -0.031 0];
                    Critic = [132.85	 34.94	 93.10  0.292  0.045];
                case 'Carbon dioxide'
                    Atoms = [1,0,2,0];
                    Prop = [-393.51	 -394.38   5.457 1.045 0 -1.157 0];
                    Critic = [304.12	 73.74	 94.07  0.274  0.225];
                case 'Ethanol'
                    Atoms = [2,6,1,0];
                    Prop = [-234.95	 -167.73   3.518 20.001 -6.002 0 0];
                    Critic = [513.92	 61.48	167.00  0.240  0.649];
                case '1-pentene'
                    Atoms = [5,10,0,0];
                    Prop = [-21.3	   78.6	   2.691 39.753 -12.447 0 0];
                    Critic = [464.8	 35.6	298.4   0.275  0.237];
                case 'Methane'
                    Atoms = [1,4,0,0];
                    Prop = [-74.52	  -50.45   4.217 1.702 9.081 -2.164 0 0];
                    Critic = [190.56   45.99   98.60  0.286  0.012];
                case 'Pentane'
                    Atoms = [5,12,0,0];
                    Prop = [-146.76	   -8.65   2.464 45.351 -14.111 0 0];
                    Critic = [469.7	 33.7	311	    0.268  0.252];
                case 'Toluene'
                    Atoms = [7,8,0,0];
                    Prop = [50.17	  122.29  0.290 47.052 -15.716 0 0];
                    Critic = [591.75	 41.08	316	    0.264  0.264];
                case 'Hydrogen'
                    Atoms = [0,2,0,0];
                    Prop = [0	    0	   3.249 0.422 0 0.083 0];
                    Critic = [32.98	 12.93	 64.2   0.303 -0.217];
                case 'Water'
                    Atoms = [0,2,1,0];
                    Prop = [-241.81	 -228.42   3.470 1.450 0 0.121 0];
                    Critic = [647.14  220.64	 55.95  0.229  0.344];
                case 'Methanol'
                    Atoms = [1,4,1,0];
                    Prop = [-200.7  -162.0 2.211 12.216 -3.450 0 0];
                    Critic = [512.64   80.97  118.00  0.224  0.565];
                case 'Ammonia'
                    Atoms = [0,3,0,1];
                    Prop = [-45.94	  -16.41  3.578 3.020 0 -0.186 0];
                    Critic = [405.4	113.53	 72.47  0.255  0.257];
                case 'Acetaldehyde'
                    Atoms = [2,4,1,0];
                    Prop = [-166.19	 -128.86  1.693 17.978 -6.158 0 0];
                    Critic = [461.2	57.0	155.5	0.300	0.310];
                case '1-butene'
                    Atoms = [4,8,0,0];
                    Prop = [-0.54	  70.37	   1.967 31.630 -9.873 0 0];
                    Critic = [419.5	40.20	240.80	0.278	0.194];
                case 'Butane'
                    Atoms = [4,10,0,0];
                    Prop = [-125.79	-16.57	   1.935 36.915 -11.402 0 0];
                    Critic = [0	0	0	0	0];
                case 'Nitrogen'
                    Atoms = [0,0,0,2];
                    Prop = [0	    0	   3.280 0.593 0 0.040 0];
                    Critic = [126.2	 33.98	 90.1   0.289  0.037];
                case 'Formaldehyde'
                    Atoms = [1,2,1,0];
                    Prop = [-116 -110 2.264 7.022 -1.877 0 0];
                    Critic = [408.0 65.90 115 0.223 0.282];
                otherwise
                    error(strcat('Nonexistent compoud ', Name, '. Add it to Compound.m AssignAtoms method'));
            end
            HoT = Prop(1)*(10^3);
            GoT = Prop(2)*(10^3);
            Cp = [Prop(3), Prop(4)*(10^-3), Prop(5)*(10^-6), Prop(6)*(10^5), Prop(7)*(10^-9)];
            Tc = Critic(1);
            Pc = Critic(2);
            Vc = Critic(3);
            Zc = Critic(4);
            Omega = Critic(5);
        end
        
    end
    
end

