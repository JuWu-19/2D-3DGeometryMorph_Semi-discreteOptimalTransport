function area=compute_surface_3d(points)
v1=points(1,:);
if size(points,1)==3
v = null(bsxfun(@minus, v1, points(2:end,:)))';
if size(v,1)>1
    area=0;
    return
end

else
ct=0;
while true
ct=ct+1;
ind_p=randsample([1:size(points,1)],3);
sub_p=points(ind_p,:);
v = null(bsxfun(@minus, sub_p(1,:), sub_p(2:end,:)))';
if size(v,1)==1||ct>5
    break;
end
end

end

if size(v,1)>1
area=0;
return
end
dotz = [0,0,1]*v';
doty = [0,1,0]*v';
dotx = [1,0,0]*v';
if dotz~=0
    if dotz < 0
    v=-1*v;
    dotz=-dotz;
    end
    
    anglez=dotz/norm(v);
    
    [~,proj_facet]=convhulln(points(:,1:2));
    
    area=proj_facet/anglez;

elseif dotx~=0
    if dotx < 0
    v=-1*v;
    dotx=-dotx;
    end
    
    anglex=dotx/norm(v);
    
    [~,proj_facet]=convhulln(points(:,2:3));
    
    area=proj_facet/anglex;

else

    if doty < 0
    v=-1*v;
    doty=-doty;
    end
    
    angley=doty/norm(v);
    
    [~,proj_facet]=convhulln(points(:,[1,3]));
    
    area=proj_facet/angley;



end

end