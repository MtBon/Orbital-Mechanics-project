function [Dvmin,row_min,col_min,at_min,et_min,ToF_min,Opt_Departure,Opt_Arrive,r1_min,r2_min,v1_min,v2_min] = FoundMin(Dvtot,ToF,...
                                                                                                   a_t,e_t,r_dep,r_arr,v_arr,v_dep,Departure, Arrive)
        
        %Departure is the linspace of the departure window;
        %arrive is the same but for the arrive;
        % c are the columns where arrives are positioned;
        % r are rows where departure are positioned;
        %The optimal datas are given in Date form;

[m,index_arr] = min(Dvtot);
[Dvmin,index_dep] = min(m);
col_min = index_dep;
row_min = (index_arr(index_dep));

ToF_min = ToF(row_min,col_min);
at_min = a_t(row_min,col_min);
et_min = e_t(row_min,col_min);

r1_min = r_dep(row_min,col_min);
r2_min = r_arr(row_min,col_min);
v1_min = v_dep(row_min,col_min);
v2_min = v_arr(row_min,col_min);

opt_dep = Departure(col_min);
opt_arr = Arrive(row_min);

Opt_Departure = mjd20002date(opt_dep);
Opt_Arrive = mjd20002date(opt_arr);


