function approximatePi = computePiPlotError(error)
    p = 0;
    sum = 0;
    index = 0;
    pApproximated = [];
    indexedError = [];
    indexedActualError = [];
    syms k;
    f = (-1)^k/(2*k + 1);
    while (abs(pi - p) > error)
        sum = sum + (-1)^index * 1/(2*index + 1);
        p = 4*sum;
        index = index + 1;
        pApproximated(index) = p;
        err = abs(pi - p);
        indexedError(index) = err;
        indexedActualError(index) = symsum(f, k, index, 10^7);
    end
    approximatePi = p;
    loglog(1:index,indexedError,'--k')
    xlabel('Iteration Index n');
    ylabel('Error');
    title('Error in Pi Approximation');
    hold on
    plot(1:index, indexedActualError, '-k')
    hold off
    legend('Error','Actual Error')
end


format long
approx_pi = computePiPlotError(5e-3)
pi - approx_pi