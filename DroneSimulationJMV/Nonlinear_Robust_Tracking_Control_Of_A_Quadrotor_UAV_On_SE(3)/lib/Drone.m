classdef Drone < handle

%% MEMBERS    
    properties
        g         % aceleracion gravitacional
        t         % tiempo de simulacion
        dt        % velocidad
        timeF     % tiempo final de la simulacion
        
      %% Parametros
        m         % masa     
        d         % Distancia desde el centro de masa al centro en el plano
        J         % Matriz de inercia con respecto al body-fixed frame
        
        %% Condiciones inciales
        s         % Variables de stado
                  % [x, v, R, omega]'
                  % [[x, y, z],[dotx, doty, dotz],[p, q, r],[phi,theta,psi]]'
        ds        % El cambio respecto a las variables de estado
                  % [dotx, dotv, dotR, dotOmega]


        x         % Vector de posicion [x, y, z]'
        v         % Vector de velocidad [dotx, doty, dotz]'
        omega     % Angulos de euler [phi, theta, psi]'
        R         % SO(3) [0, -x_3, x_2
        %                 x_3,  0, -x_1
        %                -x_2, x_1,  0]
        dotR
        dotOmega

        
      %% Entradas de control
        u         % Vector de entradas de contro [f, M1, M2, M3]'
        f         % Empuje total
        M         % Vector de momento total en body-fixed frame
                  % [M1, M2, M3]'

      
    end
    
    
    properties
        %% Parametros de control
        k_x       % Constante positiva para la posicion
        k_v       % Constante positiva para la velocidad
        k_R       % Constante positiva para matriz de rotacion
        k_Omega   % Constante positiva para vector de angulos
        c_1       % Constante positiva
        c_2       % Constante positiva
        epsilon_x % constante positiva, error de seguimiento para la posicion
        epsilon_R % Constante positiva, error de seguimiento matriz de rotacion
    end

    %% Position  Controlled Flight Mode

    properties

        position_err    % Error de posicion
        position_des    % Posicion deseada
        dotposition_des % Derivada de la posicion deseada


        R_c             %Attitude calculada
        hatOmega_c      % Velocidad angular calculada
        omega_c

        velocity_err    % Error de la velocidad
        velocity_des    % Velocidad deseada

        attitude_err    % Error de attitude
        omega_err       % Error de velocidad angular

    end

    properties
        R_d       % Acttitude deseada
        omega_d   % Velocidad angular deseada
    end
    
    
%% METHODS
    methods
    %% CONSTRUCTOR
        function obj = Drone(params, initStates, initInputs, controllerParam, simTime)

            obj.g = -9.81;
            obj.t = 0.0;
            obj.dt = 0.01;
            obj.timeF = simTime;

            %% Parametros del sistema
            % @mass = masa total del drone
            % @distance = distancia del centro de masa al centro del plano
            % @Jxx = Matriz de inercia en x
            % @Jyy = Matriz de incercia en y
            % @Jzz = Matriz de incercia en z
            % @torque = Torque generado
            
            obj.m = params('mass');
            obj.d = params('distance');
            obj.J = [params('Jxx'),0,0 ; 0,params('Jyy'),0; 0,0,params('Jzz')];
            
            %% Condiciones inciales
            % @s = Variables de estado
            % @x = Vector de posicion [x, y, z]'
            % @v = Vector de velocidad [dotx, doty, dotz]'
            % @omega = Angulos de euler [phi, theta, psi]'
            % @R = % SO(3) [0, -x_3, x_2
            %                 x_3,  0, -x_1
            %                -x_2, x_1,  0]
            
            
            obj.x = initStates('position');
            obj.v = initStates('velocity');
            obj.omega = initStates('angularV');
            obj.dotOmega = zeros(3,1);
            
            obj.R = initStates('SO(3)');
            obj.dotR = zeros(3);

            obj.s = [obj.x; obj.v; obj.omega];
            obj.ds = zeros(9, 1);
            
            obj.u = initInputs;
            obj.f = obj.u(1);
            obj.M = obj.u(2:4);
            
            %% Parametros de control
            obj.k_x = controllerParam('k_x');
            obj.k_v = controllerParam('k_v');
            obj.k_R = controllerParam('k_R');
            obj.k_Omega = controllerParam('k_Omega');
            obj.c_1 = controllerParam('c_1');
            obj.c_2 = controllerParam('c_2');
            obj.epsilon_x = controllerParam('epsilon_x');
            obj.epsilon_R = controllerParam('epsilon_R');  
        end
        
    %% Metodo para regresar las variables de estado
        function state = GetState(obj)
            state.s = obj.s;
            state.R = obj.R;
        end
        
    %% Ecuaciones dinamicas del sistema
        function obj = EvalEOM(obj)
            e_3 = [0, 0, 1]';
            
            % Translational Motions
            obj.ds(1:3) = obj.v;
            obj.ds(4:6) = 1 / obj.m * ([0; 0; obj.m * obj.g] - obj.f *obj.R * e_3);
            
            % Rotational Motions
            obj.dotR = obj.R*R2SO(obj.omega);

            obj.dotOmega = (obj.J)\ (obj.M - cross(obj.omega,obj.J*obj.omega));

        end

    %% Siguiente spacio de estados
        function obj = UpdateState(obj)
            obj.t = obj.t + obj.dt;
            
            % Find(update) the next state of obj.X
            obj.EvalEOM();
            obj.s = obj.s + obj.ds.*obj.dt;
            obj.R = obj.R + obj.dotR.*obj.dt;
            
            obj.x = obj.s(1:3);
            obj.v = obj.s(4:6);
            obj.omega = obj.s(7:9);
            obj.R = obj.R;
        end
        
    %% Controlador **** Position controlled flight mode ******************
    function obj = PositionCtrlFM(obj,refPosition)

    
        e_3 = [0,0,1]';
        b_1c = [0,1,0]';


% **************** Error de posicion *************************
        obj.position_des = refPosition;
        obj.position_err = [obj.x(1)- obj.position_des(1);
            obj.x(2) - obj.position_des(2);
            obj.x(3) - obj.position_des(3)];

% *************** Eror de velocidad *************************
        obj.dotposition_des = gradient(obj.position_des);
        obj.velocity_err =  obj.v - obj.dotposition_des;


        b_3c = - (-obj.k_x*obj.position_err- obj.k_v * obj.velocity_err - [0;0;obj.m*obj.g] + obj.m * gradient(obj.dotposition_des))/ ...
            norm(-obj.k_x*obj.position_err- obj.k_v * obj.velocity_err - [0;0;obj.m*obj.g] + obj.m * gradient(obj.dotposition_des));


        obj.R_c =[b_1c,cross(b_3c,b_1c),b_3c];
        obj.hatOmega_c = obj.R_c'*gradient(obj.R_c);

% *************** Error de attitude ************************
obj.attitude_err = (1/2)*SO2R(obj.R_c*obj.R - obj.R'*obj.R_c);

% *************** Error de velocity ************************
obj.omega_c = SO2R(obj.hatOmega_c);
obj.omega_err = obj.omega - obj.R'*obj.R_c*obj.omega_c;

aux = obj.k_x * obj.attitude_err + obj.k_v*obj.velocity_err + [0;0;obj.m*obj.g] - obj.m*gradient(obj.dotposition_des);
obj.f = dot(aux, obj.R*e_3);

obj.M = -obj.k_R*obj.attitude_err - obj.k_Omega*obj.omega_err + cross(obj.omega,obj.J*obj.omega) ...
    -obj.J*(R2SO(obj.omega)*obj.R'*obj.R_c*obj.omega_c - obj.R'*obj.R_c*gradient(obj.omega_c));








% **************** Pruebas para el lazo abierto *************************
%             obj.u(1)=-obj.m*obj.g;       
%             obj.u(2)=0.0;
%             obj.u(3)=0.0;
%             obj.u(4)=0.0;
% 
%             obj.f = obj.u(1);
%             obj.M = obj.u(2:4);
			
        end
    end
end
