function [time_serie] = Murphy2007 (maxE,maxt,observe,setvalues)
% Made by GILLESPIE_MAKER


nvar = 5;
varindex.('A') = 3;
varindex.('M') = 4;
varindex.('N') = 1;
varindex.('P') = 5;
varindex.('R') = 2;
vars = zeros(1,nvar);
vars(3) = 0; % A
vars(4) = 0; % M
vars(1) = 1; % N
vars(5) = 0; % P
vars(2) = 0; % R


nparam = 9;
paramindex.('a') = 7;
paramindex.('dmu') = 2;
paramindex.('dpii') = 4;
paramindex.('galpha') = 5;
paramindex.('grho') = 9;
paramindex.('lambda') = 8;
paramindex.('m') = 1;
paramindex.('p') = 3;
paramindex.('r') = 6;
params = zeros(1,nparam);
params(7) = 1.598000e-01; % a
params(2) = 6.931000e-01; % dmu
params(4) = 3.900000e-03; % dpii
params(5) = 2.500000e-01; % galpha
params(9) = 7.100000e-03; % grho
params(8) = 1.495000e-01; % lambda
params(1) = 10; % m
params(3) = 1; % p
params(6) = 5.370000e-02; % r


if (not(isempty(setvalues)))
    fields = fieldnames(setvalues);
    for i = 1:length(fields)
        field = fields{i};
        if (isfield(varindex,field))
            vars(varindex.(field)) = setvalues.(field);
        elseif (isfield(paramindex,field))
            params(paramindex.(field)) = setvalues.(field);
        end
    end
end

obsvars = [];
obscount = 0;
for i = 1:length(observe)
  field = observe{i};
  obscount = obscount + 1;
  obsvars(obscount) = varindex.(field);
end

time_serie = zeros(maxE,obscount+1);
nrxn = 9;
props = zeros(1,nrxn);


time = 0;
event = 0;

while((time <= maxt) && (event <= maxE))
    props(1) = params(7)*vars(1); % Act_f: N -> A
    props(2) = params(5)*vars(3); % Act_r: A -> N
    props(8) = params(2)*vars(4); % M_degrad: M -> empty
    props(9) = params(4)*vars(5); % P_degrad: P -> empty
    props(3) = params(6)*vars(1); % Rep_f: N -> R
    props(4) = params(9)*vars(2); % Rep_r: R -> N
    props(5) = params(1)*vars(3); % Transc_A: A -> A + R
    props(6) = params(8)*vars(2); % Transc_R: R -> R + M
    props(7) = params(3)*vars(4); % Transl: M -> M + P
    propsum = sum(props);
    r1 = propsum*rand(1);
    localsum = 0;
    for j = 1:nrxn
        localsum = localsum + props(j);
        if(localsum >= r1)
            if (j == 0)
            elseif (j == 1) % Act_f: N -> A
               vars(1) = vars(1)-1 ; % N-1 
               vars(3) = vars(3)+1; % A+1
               break;
            elseif (j == 2) % Act_r: A -> N
               vars(3) = vars(3)-1 ; % A-1 
               vars(1) = vars(1)+1; % N+1
               break;
            elseif (j == 8) % M_degrad: M -> empty
               vars(4) = vars(4)-1; % M-1
               break;
            elseif (j == 9) % P_degrad: P -> empty
               vars(5) = vars(5)-1; % P-1
               break;
            elseif (j == 3) % Rep_f: N -> R
               vars(1) = vars(1)-1 ; % N-1 
               vars(2) = vars(2)+1; % R+1
               break;
            elseif (j == 4) % Rep_r: R -> N
               vars(2) = vars(2)-1 ; % R-1 
               vars(1) = vars(1)+1; % N+1
               break;
            elseif (j == 5) % Transc_A: A -> A + R
               vars(4) = vars(4)+1; % M+1
               break;
            elseif (j == 6) % Transc_R: R -> R + M
               vars(4) = vars(4)+1; % M+1
               break;
            elseif (j == 7) % Transl: M -> M + P
               vars(5) = vars(5)+1; % P+1
               break;
            end
        end
    end
    time = time + (log(1/rand(1))/propsum);
    event = event + 1;
    time_serie(event,:) = [time,vars(obsvars)];
end
time_serie = time_serie(1:event,:);
end
