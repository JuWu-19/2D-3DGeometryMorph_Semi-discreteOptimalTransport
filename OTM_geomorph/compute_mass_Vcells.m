function mass=compute_mass_Vcells(PD)

num=size(PD,1);
mass=0;
for k=1:1:num
P1=polytope(PD{k});
mass=mass+volume(P1);
end

end