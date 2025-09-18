% This function draw the benchmark functions

function func_plot_cec17(func_name)

%[lb,ub,dim,fobj]=Get_Functions_details_WestPSO(func_name);
[lb,ub,dim,fobj]=CEC2017(func_name);

x=lb:2:ub;y=x;

L=length(x);
f=[];


for i=1:L
    for j=1:L
        if(length(func_name)==2)
            f(i,j)=fobj([x(i),y(j)]);
        end
        if(length(func_name)==3)
            f(i,j)=fobj([x(i),y(j),0,0,0,0,0,0,0,0]);
        end
    end
end

surfc(x,y,f,'LineStyle','none');

end

