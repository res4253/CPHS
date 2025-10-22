classdef t1dms_InputCL
    %
    properties
        input, inputArray, t_end
    end
    methods
        function self = t1dms_InputCL(t_end,input)
            self.t_end = t_end;
            self.input = input;
            self.inputArray = self.makeArray;
        end
        function input = getDeterminedInput(self,t)
            inputArray = self.inputArray;
            input = 0;
            for i = 1:size(inputArray,1)
                if t >=  inputArray(i,1) && t < inputArray(i,1) + inputArray(i,2)
                    input = inputArray(i,3) / inputArray(i,2);
                    break
                end
            end
        end

        function u = getInput(self,t)
            u = self.getDeterminedInput(t);
        end
        
        function  array = makeArray(self)
            input = self.input;
            array = [];
            sum = 0;
            count = 0;
            start_index = -1;
            for i=1:length(input)
                if input(i)==0
                    if sum > 0
                        array(end+1,:) = [start_index, count, sum];
                        sum = 0;
                        start_index = -1;
                        count = 0;
                    end
                else
                    sum = sum + input(i);
                    count = count + 1;
                    if start_index == -1
                        start_index = i-1;
                    end
                end
            end
            if sum > 0
                array(end+1,:) = [start_index, count, sum];
            end
        end
    end
end