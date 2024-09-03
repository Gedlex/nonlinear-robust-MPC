function Theta_int = intersection(Theta_1, Theta_2)
    % Intersection of two sets Theta_1 and Theta_2 in 1D
    Theta_int(1) = max([Theta_1(1) Theta_2(1)]);
    Theta_int(2) = min([Theta_1(2) Theta_2(2)]);
end