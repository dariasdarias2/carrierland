function baselineBody = baseline_creation(Nbaselines,baselineLength)
    if Nbaselines<2
        error('Not enough number of antennas!')
    end
    baselineBody = randn(3,Nbaselines);
    baselineBody = baselineLength * baselineBody./vecnorm(baselineBody, 2, 1);
    dot_products = baselineBody' * baselineBody;
    count = 0;
    while any( abs(dot_products(~eye(Nbaselines))) > 0.9 )
        baselineBody = randn(Nbaselines);
        baselineBody = baselineLength * baselineBody./vecnorm(baselineBody, 2, 1);
        dot_products = baselineBody' * baselineBody;
        count = count + 1;
        if count > 3
            break
        end
    end
end
