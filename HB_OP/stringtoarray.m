syms w cr1 cr3 ci3 ;%cr5 ci5 z;
sym_var=[w cr1 cr3 ci3];% cr5 ci5 z];

fid=fopen('N5.txt');
allData = textscan(fid,'%s','Delimiter','\n','BufSize',100000);  % bigger buffer needed
eqn={};
for j=1:length(allData{1})
    eqn{j}=char(allData{1}(j));
end
for j=1:length(sym_var),

    eqn=strrep(eqn,char(sym_var(j)),strcat('a[',int2str(j),']'));
    eqn=strrep(eqn,' ','');
end
fid=fopen('N5out.txt','w');
pre_idx=1;
for j=1:length(eqn),  % index start from pre_idx
    % add strings{i}='';
    eqn{j}=strcat(strcat('strings{',int2str(j+pre_idx-1),'}='''), char(eqn{j}), ''';\n');
    fprintf(fid,eqn{j});
end

