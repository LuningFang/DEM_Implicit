function tot_delta = getTotalPenetration(Contacts)

NC = length(Contacts);
tot_delta = 0;
for i = 1:NC
    tot_delta = tot_delta + Contacts(i).delta;
end