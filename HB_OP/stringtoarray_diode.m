syms w;
P0=sym('cr%d_0',[nsp 1]);
Pr=sym('cr%d_%d',[nsp n]);  % unknown variables real part      cr3_2: second harmonic for 3rd variable
Pi=sym('ci%d_%d',[nsp n]);  % unknown variables imag part
i_var=reshape(Pi,1,[]);%i_var(1)=[];
sym_var=[reshape(P0,1,[]), reshape(Pr,1,[]), i_var];

fid=fopen('HB_OP\N5.txt');
allData = textscan(fid,'%s','Delimiter','\n','BufSize',100000);  % bigger buffer needed
eqn={};
for j=1:length(allData{1})
    eqn{j}=char(allData{1}(j));
end
for j=1:length(sym_var),
 
    eqn=strrep(eqn,char(sym_var(j)),strcat('a[',int2str(j),']'));
    eqn=strrep(eqn,' ','');
end
fid=fopen('HB_OP\N5out.txt','w');
pre_idx=1;
for j=1:length(eqn),  % index start from pre_idx
    % add strings{i}='';
    eqn{j}=strcat(strcat('strings{',int2str(j+pre_idx-1),'}='''), char(eqn{j}), ''';\n');
    fprintf(fid,eqn{j});
end

addpath('C:\matlab\R2012a\toolbox\suiteSparse\SPQR\maTLAB\')
clear strings;