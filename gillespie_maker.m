function [data] = gillespie_maker (name)
% GILLESPIE_MAKER takes a tab delimited file (NAME.mgi) specifying the
% parameter, variables, and reactions. The syntax is as follow
% % -----------------------
% VAR    A   1
% PARAM  kf  1
% PARAM  kr  1
% REACTION	test_f
% PROP-FN	kf*A
% CHANGE-S	A-1 & B+1	
% COMMENTS	A -> B
% REACTION	test_r
% PROP-FN	kr*B
% CHANGE-S	A+1 & B-1	
% COMMENTS	B -> A
% % -----------------------
% 
% The program will create a file NAME.m that is the gillespie (SSA)
% implementation of the system described in NAME.mgi. The new
% function will take up to 3 mandatory arguments, plus one
% optional.
% MAXE := maximal number of events
% MAXT := maximal time
% OBSERVE := a cell array where each value is a string with the
% simbol of the variable you want to observe.
% SETVALUES := a structure where each field is either the string of
% a parameter of a variable. The values for each field will
% overwrite the corresponding parameter, or inital state for
% parameters and variable respectively.
% Finally this function will return a matrix where the first column
% is time, and the following are all the variables observed in the
% same order as OBSERVE.
% 
% *********************
% 
% KNOWN-BUGS: The parser is had-hock. Do not use symbol names that
% can be matched as substring of one another. The program would mix
% one and the other
% TODO: 
%  - Add more documentation
%  - Parse symbols by finding word limits
%  - ??
fid = fopen([char(name),'.mgi']);
data = textscan(fid,'%s\t%s\t%s','Delimiter','\t');
fclose(fid);

nrows = length(data{1});
nvar = 0;
nparam = 0;
nrxn = 0;
for i = 1:nrows
    col1 = data{1}(i);
    col2 = data{2}(i);
    %col3 = data{3}(i);
    if (strcmp(col1(1),'%'))
        
    elseif(strcmp(col1,'VAR'))
        nvar = nvar + 1;
        varindex.(char(col2)) = nvar;
        varinit.(char(col2)) = str2num(char(data{3}(i)));
    elseif(strcmp(col1,'PARAM'))
        nparam = nparam + 1;
        paramindex.(char(col2)) = nparam;
        paraminit.(char(col2)) = str2num(char(data{3}(i)));
    elseif(strcmp(col1,'REACTION'))
        nrxn = nrxn + 1;
        rxnindex.(char(col2)) = nrxn;
        crxn = char(col2);
    elseif(strcmp(col1,'PROP-FN'))
        rxnprop.(crxn) = char(col2);
    elseif(strcmp(col1,'CHANGE-S'))
        rxnchange.(crxn) = char(col2);
    elseif(strcmp(col1,'COMMENTS'))
        rxncomm.(crxn) = char(col2);
    else
    end
end

stream = fopen([name,'.m'],'w');
fprintf(stream,'function [time_serie] = %s (maxE,maxt,observe,setvalues)\n',name);
fprintf(stream,'%% Made by GILLESPIE_MAKER\n');
fprintf(stream,'\n\n');

vars = sort(fieldnames(varindex));
fprintf(stream,'nvar = %d;\n',nvar);
for i = 1:nvar
    var = vars{i};
    fprintf(stream,'varindex.(''%s'') = %d;\n',var,varindex.(var));
end
fprintf(stream,'vars = zeros(1,nvar);\n');
for i = 1:nvar
    var = vars{i};
    fprintf(stream,'vars(%d) = %d; %% %s\n',varindex.(var),varinit.(var),var);
end
fprintf(stream,'\n\n');

params = sort(fieldnames(paramindex));
fprintf(stream,'nparam = %d;\n',nparam);
for i = 1:nparam
    param = params{i};
    fprintf(stream,'paramindex.(''%s'') = %d;\n',param,paramindex.(param));
end
fprintf(stream,'params = zeros(1,nparam);\n');
for i = 1:nparam
    param = params{i};
    fprintf(stream,'params(%d) = %d; %% %s\n',paramindex.(param),paraminit.(param),param);
end
fprintf(stream,'\n\n');

fprintf(stream,'if (not(isempty(setvalues)))\n');
fprintf(stream,'    fields = fieldnames(setvalues);\n');
fprintf(stream,'    for i = 1:length(fields)\n');
fprintf(stream,'        field = fields{i};\n');
fprintf(stream,'        if (isfield(varindex,field))\n');
fprintf(stream,'            vars(varindex.(field)) = setvalues.(field);\n');
fprintf(stream,'        elseif (isfield(paramindex,field))\n');
fprintf(stream,'            params(paramindex.(field)) = setvalues.(field);\n');
fprintf(stream,'        end\n');
fprintf(stream,'    end\n');
fprintf(stream,'end\n\n');

fprintf(stream,'obsvars = [];\n');
fprintf(stream,'obscount = 0;\n');
fprintf(stream,'for i = 1:length(observe)\n');
fprintf(stream,'  field = observe{i};\n');
fprintf(stream,'  obscount = obscount + 1;\n');
fprintf(stream,'  obsvars(obscount) = varindex.(field);\n');
fprintf(stream,'end\n\n');

fprintf(stream,'time_serie = zeros(maxE,obscount+1);\n');

rxns = sort(fieldnames(rxnindex));
fprintf(stream,'nrxn = %d;\n',nrxn);
fprintf(stream,'props = zeros(1,nrxn);\n');
fprintf(stream,'\n\n');

fprintf(stream,'time = 0;\n');
fprintf(stream,'event = 0;\n\n');
fprintf(stream,'while((time <= maxt) && (event <= maxE))\n');
for i = 1:nrxn
    rxn = rxns{i};
    prop_fn = replace_into_string(rxnprop.(rxn));
    fprintf(stream,'    props(%d) = %s; %% %s: %s\n',rxnindex.(rxn),prop_fn,rxn,rxncomm.(rxn));
end
fprintf(stream,'    propsum = sum(props);\n');
fprintf(stream,'    r1 = propsum*rand(1);\n');
fprintf(stream,'    localsum = 0;\n');
fprintf(stream,'    for j = 1:nrxn\n');
fprintf(stream,'        localsum = localsum + props(j);\n');
fprintf(stream,'        if(localsum >= r1)\n');
fprintf(stream,'            if (j == 0)\n');
for i = 1:nrxn
    rxn = rxns{i};
    fprintf(stream,'            elseif (j == %d) %% %s: %s\n',rxnindex.(rxn),rxn,rxncomm.(rxn));
    updates = textscan(rxnchange.(rxn),'%s','Delimiter','&');
    updates = updates{1};
    for j = 1:length(updates)
        update = updates{j};
        [right,left] = replace_into_string(update);
        fprintf(stream,'               %s = %s; %% %s\n',left,right,update);
    end
    fprintf(stream,'               break;\n');
end

fprintf(stream,'            end\n');
fprintf(stream,'        end\n');
fprintf(stream,'    end\n');
fprintf(stream,'    time = time + (log(1/rand(1))/propsum);\n');
fprintf(stream,'    event = event + 1;\n');
fprintf(stream,'    time_serie(event,:) = [time,vars(obsvars)];\n');
fprintf(stream,'end\n');
fprintf(stream,'time_serie = time_serie(1:event,:);\n');
fprintf(stream,'end\n');
fclose(stream);

function [result,last] = replace_into_string (s) % see bugs
    result = '';
    temps = s;
    while(~isempty(temps))
        next_symb = '';
        next_symb_val = '';
        for q = 1:length(temps)
            symb = temps(1:q);
            if (isfield(varindex,symb))
                next_symb = char(symb);
                next_symb_val = ['vars(',num2str(varindex.(symb)),')'];
                break;
            elseif (isfield(paramindex,symb))
                next_symb = char(symb);
                next_symb_val = ['params(',num2str(paramindex.(symb)),')'];
                break;
            end
        end
        %next_symb
        %next_symb_val
        if (isempty(next_symb))
            result = [result,temps(1)];
            if (length(temps) > 1)
                temps = temps(2:end);
            else
                temps = '';
            end
        else
            result = [result,next_symb_val];
            last = next_symb_val;
            if (length(temps) > length(next_symb))
                temps = temps(length(next_symb)+1:end);
            else
                temps = '';
            end
        end
    end
end

end