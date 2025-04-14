function L=  BC_periodic_IMEX_L(L,D)


%periodic boundary
L(1,:)=0;
L(1,1)=1;
L(1,end)=-1;% y(1)=y(end);

%dy=D*y, dy(1)=dy(end)
L(end,:)=D(1,:)-D(end,:);


% L(2,:)=1;
end