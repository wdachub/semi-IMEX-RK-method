function dy = myodeIM(t,y)
dy=spdiags(-y,0,length(y),length(y))
end