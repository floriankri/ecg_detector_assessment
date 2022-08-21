function P = gen_PPGpulse(t, PPGpulseType)
    
    % 1 x Lognromal & 2 x Gauss functions
    g0 = ( PPGpulseType.k0 * lognpdf((t - PPGpulseType.tau0), PPGpulseType.mu0, PPGpulseType.sigma0) );
    g1 = ( PPGpulseType.k1 * gaussmf(t,                      [PPGpulseType.sigma1 PPGpulseType.mu1]) );
    g2 = ( PPGpulseType.k2 * gaussmf(t,                      [PPGpulseType.sigma2 PPGpulseType.mu2]) );
    
    % Sum of first 1 x Lognormal & 2 x Gauss functions
    P = g0 + g1 + g2;
end

    