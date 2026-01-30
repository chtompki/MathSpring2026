function approximatePi = computePiPlotError(error)
    p = 0;
    sum = 0;
    index = 0;
    pApproximated = [];
    indexedError = [];
    while (abs(pi - p) > error)
        sum = sum + (-1)^index * 1/(2*index + 1);
        p = 4*sum;
        index = index + 1;
        pApproximated(index) = p;
        err = abs(pi - p);
        indexedError(index) = err;
    end
    approximatePi = p;
    loglog(1:index,indexedError,'--k')
    xlabel('Iteration Index n');
    ylabel('Error');
    title('Error in Pi Approximation');
end


format long
approx_pi = computePiPlotError(5e-3)
pi - approx_pi