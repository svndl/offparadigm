function template = makeElement(baseSize, baseLum, pedestalSize, pedestalLum, probeSize, probeLum)

    % do a 1/4 template and mirror-flip twice
    probe_first_half = round(probeSize*.5);
    probe_second_half = probeSize - probe_first_half;
    
    pedestal_first_half =  round(pedestalSize*.5);
    pedestal_second_half = pedestalSize - pedestal_first_half;
    
    base_first_half = round(baseSize*.5);
    base_second_half = baseSize - base_first_half;
    
    % first q
    
    template_q1 = baseLum*ones(base_first_half, base_first_half);
    template_q1(1:pedestal_first_half, 1:pedestal_first_half) = pedestalLum;
    template_q1(1:probe_first_half, 1:probe_first_half) = probeLum;

    % second q
    template_q2 = baseLum*ones(base_first_half, base_second_half); 
    template_q2(1:pedestal_first_half, 1:pedestal_second_half) = pedestalLum;
    if (probe_second_half)
        template_q2(1:probe_second_half, 1:probe_second_half) = probeLum;
    end
    
    % third q
    template_q3 = baseLum*ones(base_second_half, base_first_half);
    template_q3(1:pedestal_second_half, 1:pedestal_first_half) = pedestalLum;
    
    if (probe_second_half)
        template_q3(1:probe_second_half, 1:probe_first_half) = probeLum;
    end
    % q4
    template_q4 = baseLum*ones(base_second_half, base_second_half);
    template_q4(1:pedestal_second_half, 1:pedestal_second_half) = pedestalLum;

    if (probe_second_half)
        template_q4(1:probe_second_half, 1:probe_second_half) = probeLum;
    end
    template = [template_q1(end:-1:1, end:-1:1), template_q2(end:-1:1, :); template_q3(:, end:-1:1), template_q4];
end