function filtered_vec = filter_spacing(vec, min_dist)
    vec = sort(vec); % Ensure the vector is sorted
    filtered_vec = vec(1); % Keep the first element
    
    for i = 2:length(vec)
        if vec(i) - filtered_vec(end) >= min_dist
            filtered_vec = [filtered_vec, vec(i)]; % Append the element if it meets the spacing criterion
        end
    end
end
