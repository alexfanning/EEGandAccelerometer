

function dataOut = taskList(parameters,dataOut,dataStruct,type)

a = 1;
% Loop through each task to see if it exists in subjects data
for i = 1:length(parameters.sheetNames)
    for m = 1:length(parameters.names)
        
        tf = strcmpi(parameters.sheetNames{i},parameters.names{m});
        if type == 1
            
            if tf == 1 && i == 1
                dataOut{1,1} = dataStruct.Action.F;
            end
    
            if tf == 1 && i == 2
                dataOut{1,2} = dataStruct.Posture.F;
            end
    
            if tf == 1 && i == 3
                dataOut{1,3} = dataStruct.Spiral.F;
            end
    
            if tf == 1 && i == 4
                dataOut{1,4} = dataStruct.Tapping.F;
            end
    
            if tf == 1 && i == 5
                dataOut{1,5} = dataStruct.Rest.F;
            end

        elseif type == 2

            if tf == 1 && i == 1
                dataOut.taskNums(m) = 1;
            end
    
            if tf == 1 && i == 2
                dataOut.taskNums(m) = 2;
            end
    
            if tf == 1 && i == 3
                dataOut.taskNums(m) = 3;
            end
    
            if tf == 1 && i == 4
                dataOut.taskNums(m) = 4;
            end
    
            if tf == 1 && i == 5
                dataOut.taskNums(m) = 5;
            end
        end
    end
end

if type == 2
    dataOut.taskNums(2) = 6;
end