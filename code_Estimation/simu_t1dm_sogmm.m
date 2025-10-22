function data = simu_t1dm_sogmm(Ts,t_end,options)
    % 出力データ
    % data.ts = ts; 
    % data.xs = xs'; 
    % data.sogmm = sogmm;
    % data.magni_model = magni_model;

    arguments
        Ts
        t_end
        options.sogmm = []
        options.insulin = []
        options.meal = []
    end
  
    sogmm = options.sogmm;
    insulin = options.insulin;
    meal = options.meal;

    mealInput = t1dms_InputCL(t_end,meal);
    insulinInput = t1dms_insulinCL(insulin);

    odeoption = odeset("RelTol",1e-6,"MaxStep",0.5);
    [ts,xs] = ode45(@(t,x)sogmm.dynamics(t,x,...
        [mealInput.getInput(t);(1/6)*insulinInput.getInsulinInput(t)],sogmm.p_0),...
        0:Ts:t_end,sogmm.x_0,odeoption);


    data.ts = ts'; 
    data.xs = xs'; 
    data.Gs = xs(:,1)';
    data.sogmm = sogmm;
    data.array = insulinInput;

    i=1;
    for t=0:Ts:t_end
        data.meal(i)=mealInput.getInput(t);
        i=i+1;
    end
    i=1;
    for t=0:Ts:t_end
        data.insulin(i)=(1/6)*insulinInput.getInsulinInput(t);
        i=i+1;
    end

end