classdef DsCL < t1dms_InputCL
    properties
        Ds , Ds_0  , q_0 , DsArray , GastrointestinalSubmodel
    end
    methods
        function self = DsCL(GastrointestinalSubmodel,option)
            arguments
                GastrointestinalSubmodel
                option.input
                option.t_end (1,1) {mustBeNumeric} = 99999999;
                option.array =[];
            end
            self@t1dms_InputCL(option.t_end,option.input)
            self.Ds_0 = 0; %0か45000　どっち?
            self.GastrointestinalSubmodel = GastrointestinalSubmodel;
            self.q_0 = GastrointestinalSubmodel.x_0;

            if isempty(option.array)
                self.inputArray = self.makeArray;%data[0 0 0]
            else
                self.inputArray = option.array;
            end

            self.DsArray = self.makeDsArray;
            
        end

        function DsArray = makeDsArray(self,~)
            inputArray = self.inputArray;
            q_0 = self.q_0;

            DsArray = [];
            for phaseNum = 1:size(inputArray,1)
                Ds = q_0(1) + q_0(2) + inputArray(phaseNum,3);

                if phaseNum == size(inputArray,1)
                    DsArray = [DsArray;
                               inputArray(phaseNum,1), self.t_end-inputArray(phaseNum,1),Ds];
                    break
                else
                    DsArray = [DsArray;
                              inputArray(phaseNum,1), inputArray(phaseNum+1,1)-inputArray(phaseNum,1),Ds];
                end

                t1 = inputArray(phaseNum,1);
                t2 = inputArray(phaseNum+1,1);
                options = odeset('MaxStep', 0.5, 'RelTol',1e-12);
                [t, q] = ode45(@(t,q)self.GastrointestinalSubmodel.dynamics(t,q,0,[self.getInput(t),Ds]),[t1,t2],q_0, options);
                q_0 = q(end,:);
                
            end
        end

        function input = getDsDeterminedInput(self, t) 
            inputArray = self.DsArray;
            input = self.Ds_0;
            for i = 1:size(inputArray,1)
                if t >=  inputArray(i,1) && t < inputArray(i,1) + inputArray(i,2)
                    input = inputArray(i,3);
                    break
                end
            end
        end

        function Ds = getDsInput(self, t)
            Ds = self.getDsDeterminedInput(t);
        end
    end
end
               



                

