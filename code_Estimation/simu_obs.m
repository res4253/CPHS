function data = simu_obs(Ts,t_end,t_span,mode,G,L,options)
    % 出力データ
    % data.ts = ts; 
    % data.xs = xs'; 
    % data.sogmm = sogmm;
    % data.magni_model = magni_model;

    arguments
        Ts
        t_end
        t_span
        mode
        G
        L
        options.ssogmm = []
        options.insulin = []
        options.meal = []
        options.dynamics
        options.q =[];
    end
  
    ssogmm = options.ssogmm;
    insulin = options.insulin;
    meal = options.meal;

    mealInput = t1dms_InputCL(t_end,meal);
    insulinInput = t1dms_insulinCL(insulin);
    % 
    % submodel = GastrointestinalSubmodelCL(ssogmm.p_0,options.q);
    % Ds = DsCL(submodel,'t_end',t_end,'input',meal,'array',[]);
    % 
    % dynamics = options.dynamics;

    odeoption = odeset("RelTol",1e-6,"MaxStep",0.5);


    [ts,xs] = ode45(@(t,x)obs(t,x,...
        [mealInput.getInput(t);(1/6)*insulinInput.getInsulinInput(t);get_input(t,Ts,mode);get_input(t,Ts,G)],ssogmm.p_0,L),...
        t_span,ssogmm.x_0,odeoption);


    data.ts = ts'; 
    data.xs = xs'; 
    data.Gs = xs(:,1)';
    data.ssogmm = ssogmm;
    data.array = insulinInput;

    i=1;
    for t=t_span
        data.meal(i)=mealInput.getInput(t);
        i=i+1;
    end
    i=1;
    for t=t_span
        data.insulin(i)=(1/6)*insulinInput.getInsulinInput(t);
        i=i+1;
    end

    i=1;
    for t=t_span
        data.modes(i) = get_input(t,Ts,mode);
        i=i+1;
    end

    function input = get_input(t,Ts,array)
        id = floor((t-t_span(1))/Ts) + 1;%floor その要素以下の最も近い整数に丸める
        input = array(:,id);
    end
end