function FLAG = check_zero_repeat_overlap(cord, bounds, r)

FLAG = 0;

% Validate that there are no particle repeats (i.e., complete overlap).
% If there are particle repeats, remove all but one of the repeat particles
% and send an error to the user that the simulation may be incorrect. 
dummy = unique(cord,'rows'); %Breaks 3D necessary for 2D
if size(dummy,1) ~= size(cord,1)
    FLAG = 1;
    disp('ERROR! Particles were repeated.')
    disp('Repeat particles are being removed.')
    disp('Fill fraction and other metrics may be incorrect!')
    cord = dummy;
end
[flag,loc] = ismember([0,0,0],cord,'rows');
if flag == 1
    cord(loc,:) = [];
    cord = [0,0,0;cord];
else
    FLAG = 1;
    disp('ERROR! No particle was generated in the center!')
    %keyboard;
end

% Check and make sure no particles are overlapping. 
[overlap_distance, overlap_idx] = ...
    check_distance_function_v2(cord, r, bounds(1,:), bounds(2,:));

if sum(overlap_idx) ~= 0
    disp(['ERROR! Particles created are overlaping!'])
    disp(['This may break other functions and return wrong FF metrics'])
    FLAG = 1;
end

end

