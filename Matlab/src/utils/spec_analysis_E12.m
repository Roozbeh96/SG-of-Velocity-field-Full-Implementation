function spec_analysis_E12(obj)
    %{
    More info about how to compute the spectra in:
    https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/
    article/stochastic-modal-velocity-field-in-roughwall-turbulence/
    800D97E60609F4FE4D92847CCDB1C7A0
    %}
    fprintf('Spectral Analysis E12...\n');
    uprime = obj.Gen_u_HRVF - mean(obj.Gen_u_HRVF,3);
    wprime = obj.Gen_w_HRVF - mean(obj.Gen_w_HRVF,3);
    concat_uprime = reshape(uprime,size(uprime, 1), []);
    concat_wprime = reshape(wprime,size(wprime, 1), []);
    ksmax=floor(2*pi/(obj.Gen_x_HRVF(2)-obj.Gen_x_HRVF(1)));
    
    win15del=floor(15*obj.delta/(obj.Gen_x_HRVF(2)-obj.Gen_x_HRVF(1)));
    window=hann(win15del);
    NFFT=2^(nextpow2(win15del));
    [CPSD_k,k]=cpsd(concat_uprime',concat_wprime',window,NFFT/2,win15del,ksmax); %PSD_f[m^2/s]
    Apk = trapz(k,real(CPSD_k));
    Cov = sum(uprime.*wprime,3)/(size(uprime,3)-1);
    Q = mean(Cov,2)./Apk';  %%[m^3/s^2]=[m^2/s]*[m^2/s^2]*1/[m/s]
    obj.E12=(real(CPSD_k).*Q')';
end