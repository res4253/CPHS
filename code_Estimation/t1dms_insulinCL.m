%かなり無理やり
classdef t1dms_insulinCL
    properties
        input, inputArray, ins_0,test
    end
    methods 
        function self = t1dms_insulinCL(input)
            self.input = input;
            self.ins_0 = self.input(1);
            self.inputArray = self.makeBolusArray;
        end

        function array = makeBolusArray(self)

            ins_0 = self.ins_0;
            x = self.input;
            input = x - ins_0;

            n = 0;% 後方に拡張する個数

            idxs = find(input >= 6000);% 値が6000以上の位置を取得

            for i = 1:length(idxs)
                start_idx = idxs(i);
                end_idx = min(length(input), start_idx + n);

                input(start_idx:end_idx) = input(start_idx)/(n+1);  % 同じ値で上書き
            end


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

        function insulinInput = getInsulinInput(self, t)
            determinedInsulinInput= self.getDeterminedInput(t);
            insulinInput = determinedInsulinInput + self.ins_0;
        end

    end
end




