classdef Reactor
    %Reactor
    %   Detailed explanation goes here
    
    properties
        Reactions            % Array of Reaction objects
        NReactions           % Number of reactions
        N0                   % Array with nitial number of moles of each compound
        P                    % Pressure of the system
        T0                   % Initial temperature of the system
        Tf                   % Final temperature of the system
        Keq                  % Array with the reaction equilibrium constants
    end
    
    methods
        
        function obj = Reactor(Reactions, N0, P, T0, Tf)
            %Reactor Construct an instance of this class
            
            obj.Reactions = Reactions;
            obj.NReactions = length(Reactions);
            obj.N0 = N0;
            obj.P = P;
            obj.T0 = T0;
            obj.Tf = Tf;
            obj.Keq = CalculateKeq(obj, Reactions);
        end
        
    end
    
    methods (Access = private)
        
        function Keq = CalculateKeq(obj, Reactions)
            % Por enquanto calcular só do primeiro
            rxntest = Reactions(1).Compounds
            keyboard
        end
        
    end
end

