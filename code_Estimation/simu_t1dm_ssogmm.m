function data = simu_t1dm_ssogmm(Ts,t_end,t_span,mode,options)
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

    submodel = GastrointestinalSubmodelCL(ssogmm.p_0,options.q);
    Ds = DsCL(submodel,'t_end',t_end,'input',meal,'array',[]);


    dynamics = options.dynamics;

    odeoption = odeset("RelTol",1e-6,"MaxStep",0.5);

    if dynamics == "ssogmm"
        [ts,xs] = ode45(@(t,x)ssogmm_dynamics(t,x,...
            [mealInput.getInput(t);(1/6)*insulinInput.getInsulinInput(t);get_input(t,Ts,mode)],ssogmm.p_0),...
            t_span,ssogmm.x_0,odeoption);
    elseif dynamics == "inequality"
        temp = 0:Ts:t_end;
        [ts,xs] = ode45(@(t,x)inequality_dynamics(t,x,...
            [mealInput.getInput(t);(1/6)*insulinInput.getInsulinInput(t);Ds.getDsInput(t)],ssogmm.p_0,options.q),...
            temp,ssogmm.x_0,odeoption);
    elseif dynamics == "linear"
        [ts,xs] = ode45(@(t,x)linear_dynamics(t,x,...
            [mealInput.getInput(t);(1/6)*insulinInput.getInsulinInput(t);get_input(t,Ts,mode)],ssogmm.p_0),...
            t_span,ssogmm.x_0,odeoption);
    end



    data.ts = ts'; 
    data.xs = xs'; 
    data.Gs = xs(:,1)';
    data.ssogmm = ssogmm;
    data.insulinCL = insulinInput;
    data.DsCL = Ds;

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
    for t=0:Ts:t_end
        data.Ds(i)=Ds.getDsInput(t);
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